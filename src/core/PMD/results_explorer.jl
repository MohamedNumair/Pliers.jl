#=
The Idea here is to get the results dictionary of PMD power flow or optimal power flow and then match it with eng files.
=#


"""
calc_bases_from_dict(data::Dict{String,Any}; return_dict=false) 

Calculate electrical base quantities from a data dictionary.

# Arguments
- `data::Dict{String,Any}`: A dictionary containing network data. It must include:
    - `"per_unit"`: (optional) A boolean indicating whether the values are in per-unit.
    - `"settings"`: A nested dictionary with the following keys:
        - `"vbases_default"`: A collection where the first element has a `second` field representing the base voltage in volts.
        - `"voltage_scale_factor"`: A scaling factor for voltage.
        - `"sbase_default"`: The default base apparent power in VA.
        - `"power_scale_factor"`: A scaling factor for power.
- `return_dict::Bool=false`: If `true`, returns a dictionary of base quantities; otherwise, returns a tuple.

# Returns
If `return_dict == false`, returns a tuple:
1. `is_perunit::Bool`: Indicates if the values are in per-unit.
2. `vbase_V::Float64`: The base voltage in volts (phase-to-neutral).
3. `sbase_VA::Float64`: The base apparent power in volt-amperes.
4. `Zbase_Ω::Float64`: The base impedance in ohms.
5. `Ibase_A::Float64`: The base current in amperes (phase current).
6. `vbase_ll::Float64`: The base line-to-line voltage in volts.
7. `Ibase_A_ll::Float64`: The base line current in amperes.
8. `Ibase_A_ϕ::Float64`: The base phase current in amperes.

If `return_dict == true`, returns a `Dict{String,Any}` with the following keys:
- `"is_perunit"`
- `"vbase_V"`
- `"sbase_VA"`
- `"Zbase_Ω"`
- `"Ibase_A"`
- `"vbase_ll"`
- `"Ibase_A_ll"`
- `"Ibase_A_ϕ"`

# Example
- is_perunit, vbase_V, sbase_VA, Zbase_Ω, Ibase_A, vbase_ll, Ibase_A_ll, Ibase_A_ϕ = calc_bases_from_dict(data)
- bases = calc_bases_from_dict(data; return_dict=true)
"""
function calc_bases_from_dict(data::Dict{String,Any}; return_dict = false)   
    

    is_perunit = haskey(data, "per_unit") ? data["per_unit"] : false
    vbase_V = first(data["settings"]["vbases_default"]).second*data["settings"]["voltage_scale_factor"]
    vbase_ll = vbase_V*sqrt(3)
    sbase_VA = data["settings"]["sbase_default"]*data["settings"]["power_scale_factor"]
    Zbase_Ω = vbase_ll^2/sbase_VA
    Ibase_A = sbase_VA/vbase_V

    Ibase_A_ll = sbase_VA/(sqrt(3)*vbase_ll)
    Ibase_A_ϕ = sbase_VA/(3*vbase_V)

    if return_dict
        bases = Dict()
        bases["is_perunit"] = is_perunit
        bases["vbase_V"] = vbase_V
        bases["sbase_VA"] = sbase_VA
        bases["Zbase_Ω"] = Zbase_Ω
        bases["Ibase_A"] = Ibase_A
        bases["vbase_ll"] = vbase_ll
        bases["Ibase_A_ll"] = Ibase_A_ll
        bases["Ibase_A_ϕ"] = Ibase_A_ϕ
        return bases
    end

    return is_perunit, vbase_V, sbase_VA, Zbase_Ω, Ibase_A, vbase_ll, Ibase_A_ll, Ibase_A_ϕ
end



# function to get the results dictionary

function pf_results(eng::Dict{String, Any}; max_iter = 5000,kwargs...)
    # I won't develop this function further, the strategy will be to calcualte your results using typical PMD.jl functions and then you use dictify solution and viz anything. 
    is_explicit_netrual = length(eng["conductor_ids"]) > 3 ? true :  false

    _is_eng(eng) ? math=transform_data_model(eng, kron_reduce=!is_explicit_netrual, phase_project=!is_explicit_netrual) : error("This function only supports ENGINEERING data model for the moment")
    
    add_start_vrvi!(math)
    PF = compute_mc_pf(math; explicit_neutral=is_explicit_netrual, max_iter=max_iter)

    return pf_results(PF,math, eng; kwargs...)
end



function pf_results(results::Dict{String, Any}, math::Dict{String, Any}, eng::Dict{String, Any}; detailed = false, keep_pu = false)
    
    # This is just the pretty table summary for results
    results_headings = ["Termination Status", "Iterations", "Total Time", "Build Time", "Post Time", "Solve Time", "Stationarity"]
    results_values = [results["termination_status"], results["iterations"], results["time_total"], results["time_build"], results["time_post"], results["time_solve"], results["stationarity"]]
    summary_table = DataFrame(Label = results_headings, Value = results_values)
    pretty_table(summary_table, highlighters=(_highlight_results_status))
    
    
    # Moving results into the ENG dictionary 
    eng = deepcopy(eng)
    pf_sol = results["solution"]
    is_perunit, vbase_V, sbase_VA, Zbase_Ω, Ibase_A, vbase_ll, Ibase_A_ll, Ibase_A_ϕ  = calc_bases_from_dict(pf_sol)

    for (b, bus) in pf_sol["bus"]
        bus["V"] = bus["vm"] .* exp.(im * bus["va"])
        terminals = math["bus"][b]["terminals"]

        # write a dictionary where the key is the terminal number and the value is the voltage at that terminal
        bus["V"] = Dict(string(term) => bus["V"][i] for (i, term) in enumerate(terminals))

        eng_bus_id = math["bus"][b]["name"]

        if eng_bus_id == "_virtual_bus.voltage_source.source"
            eng["virtual_bus"] = Dict()
            eng["virtual_bus"]["V"] = bus["V"]
            eng["virtual_bus"]["math_id"] = b
        else
            eng["bus"][eng_bus_id]["V"] = bus["V"]
            eng["bus"][eng_bus_id]["math_id"] = b
        end
    end

    for (br, branch) in pf_sol["branch"]
        len = Int(length(branch["cr"]) / 2)
        I_t = branch["cr"][1:len] + im * branch["ci"][1:len]
        I_f = branch["cr"][len+1:end] + im * branch["ci"][len+1:end]

        branch["I_t"] = Dict(string(term) => I_t[i] for (i, term) in enumerate(math["branch"][br]["t_connections"]))
        branch["I_f"] = Dict(string(term) => I_f[i] for (i, term) in enumerate(math["branch"][br]["f_connections"]))


        f_bus = string(math["branch"][br]["f_bus"])
        t_bus = string(math["branch"][br]["t_bus"])

        Vf = pf_sol["bus"][f_bus]["V"]
        Vt = pf_sol["bus"][t_bus]["V"]

        branch["S_f"] = Dict(string(term) => Vf[string(term)] * conj(branch["I_f"][string(term)]) for (_, term) in enumerate(math["branch"][br]["f_connections"]))
        branch["S_t"] = Dict(string(term) => Vt[string(term)] * conj(branch["I_t"][string(term)]) for (_, term) in enumerate(math["branch"][br]["t_connections"]))

        branch["S_f"] = Dict(string(term) => Vf[string(term)] * conj(branch["I_f"][string(term)]) for (_, term) in enumerate(math["branch"][br]["f_connections"]))
        branch["S_t"] = Dict(string(term) => Vt[string(term)] * conj(branch["I_t"][string(term)]) for (_, term) in enumerate(math["branch"][br]["t_connections"]))

        eng_branch_id = math["branch"][br]["name"]

        if eng_branch_id == "_virtual_branch.voltage_source.source"
            eng["virtual_branch"] = Dict()
            eng["virtual_branch"]["I_t"] = branch["I_t"]
            eng["virtual_branch"]["I_f"] = branch["I_f"]
            eng["virtual_branch"]["S_t"] = branch["S_t"]
            eng["virtual_branch"]["S_f"] = branch["S_f"]
            eng["virtual_branch"]["V_f"] = pf_sol["bus"][f_bus]["V"]
            eng["virtual_branch"]["V_t"] = pf_sol["bus"][t_bus]["V"]   
            eng["virtual_branch"]["math_id"] = br
        else
            eng["line"][eng_branch_id]["I_t"] = branch["I_t"]
            eng["line"][eng_branch_id]["I_f"] = branch["I_f"]
            eng["line"][eng_branch_id]["S_t"] = branch["S_t"]
            eng["line"][eng_branch_id]["S_f"] = branch["S_f"]
            eng["line"][eng_branch_id]["V_f"] = pf_sol["bus"][f_bus]["V"]
            eng["line"][eng_branch_id]["V_t"] = pf_sol["bus"][t_bus]["V"]
            eng["line"][eng_branch_id]["math_id"] = br
        end
    end

    for (g, gen) in pf_sol["gen"]
        Ig = gen["crg"] + im * gen["cig"]

        gen["Ig"] = Dict(string(term) => Ig[i] for (i, term) in enumerate(math["gen"][g]["connections"]))

        gen_bus = string(math["gen"][g]["gen_bus"])
        Vg = pf_sol["bus"][gen_bus]["V"]

        gen["Sg"] = Dict(string(term) => Vg[string(term)] * conj(gen["Ig"][string(term)]) for (i, term) in enumerate(math["gen"][g]["connections"]))

        gen_eng_id = math["gen"][g]["name"]

        if gen_eng_id == "_virtual_gen.voltage_source.source"
            #eng["voltage_source"]["source"] = Dict()
            eng["voltage_source"]["source"]["I_g"] = gen["Ig"]
            eng["voltage_source"]["source"]["S_g"] = gen["Sg"]
            eng["voltage_source"]["source"]["math_id"] = g
        else
            eng["voltage_source"][gen_eng_id]["I_g"] = gen["Ig"]
            eng["voltage_source"][gen_eng_id]["S_g"] = gen["Sg"]
            eng["voltage_source"][gen_eng_id]["math_id"] = g
        end
    end

    for (l, load) in pf_sol["load"]
        Il = load["crd"] + im * load["cid"]

        load["Il"] = Dict(string(term) => Il[i] for (i, term) in enumerate(math["load"][l]["connections"]))

        load_bus = string(math["load"][l]["load_bus"])
        Vl = pf_sol["bus"][load_bus]["V"]

        load["Sl"] = Dict(string(term) => Vl[string(term)] * conj(load["Il"][string(term)]) for (i, term) in enumerate(math["load"][l]["connections"]))

        load_eng_id = math["load"][l]["name"]
        eng["load"][load_eng_id]["I_l"] = load["Il"]
        eng["load"][load_eng_id]["S_l"] = load["Sl"]
        eng["load"][load_eng_id]["math_id"] = l
    end

    eng["bases"] = Dict(
        "is_perunit" => is_perunit,
        "vbase_V" => vbase_V,
        "sbase_VA" => sbase_VA,
        "Zbase_Ω" => Zbase_Ω,
        "Ibase_A" => Ibase_A,
        "vbase_ll" => vbase_ll,
        "Ibase_A_ll" => Ibase_A_ll,
        "Ibase_A_ϕ" => Ibase_A_ϕ
    )

    if detailed
        pf_results_buses(pf_sol, math; keep_pu)
    end
    
    return eng, math, results
end



function branch_viz(pf_sol::Dict{String, Any}, math, branch_id; keep_pu=true::Bool, makie_backend=WGLMakie, rounding_digits = 4::Int, fig_size= nothing::Union{Tuple{Int, Int}, Nothing}, eng=nothing::Union{Dict{String, Any},Nothing})
   
    makie_backend.activate!()
    dictify_solution!(pf_sol, math)
    is_perunit, vbase_V, sbase_VA, Zbase_Ω, Ibase_A, vbase_ll, Ibase_A_ll, Ibase_A_ϕ = calc_bases_from_dict(pf_sol)
    terms = length(math["branch"][branch_id]["f_connections"])
    fig_size = isnothing(fig_size) ? (200*2*terms, 200*terms) : fig_size
    f = Figure(size = fig_size)

    f_bus_math_id = math["branch"][branch_id]["f_bus"]
    t_bus_math_id = math["branch"][branch_id]["t_bus"]
    
    branch_c_fr = pf_sol["branch"][branch_id]["current_from"]
    branch_c_to = pf_sol["branch"][branch_id]["current_to"]
    
    branch_imp_matrix = math["branch"][branch_id]["br_r"] .+ im * math["branch"][branch_id]["br_x"]


    gl_fbus = GridLayout(f[1,1]) # Left panel for 'from' bus 
    gl_mid = GridLayout(f[1,2]) # Middle panel for branch details
    gl_tbus = GridLayout(f[1,3]) # Right panel for 'to' bus

    eng_line_id = math["branch"][branch_id]["name"]


    if isnothing(eng) || occursin("virtual", eng_line_id)
        linecode = "N/A" 
        line_len = "N/A"
        @info "if you want linecode data pass the ENGINEERING model with the eng argument"
        lc_imp = fill("", size(branch_imp_matrix))
    else
         linecode = eng["line"][eng_line_id]["linecode"]
         line_len = eng["line"][eng_line_id]["length"]
         lc_imp = eng["linecode"][linecode]["rs"] .+ im * eng["linecode"][linecode]["xs"]
         lc_imp = keep_pu ? lc_imp : lc_imp * Zbase_Ω;
    end 
    
    Label(gl_mid[0,1], "Branch: $branch_id [$(eng_line_id)] lc: [$linecode] length: [$line_len] m", fontsize=16, tellwidth=false, justification=:center,color=:blue)
    
    
    gl_matrix = GridLayout(gl_mid[1,1])
    for i in axes(branch_imp_matrix, 1)
        for j in axes(branch_imp_matrix, 2)            
                    if !ismissing(branch_imp_matrix[i, j])
                        branch_imp_matrix[i, j] = keep_pu ? branch_imp_matrix[i, j] : branch_imp_matrix[i, j] * Zbase_Ω;
                        zij = round(branch_imp_matrix[i, j], digits=rounding_digits);
                        i ==j ? _MK.Box(gl_matrix[i, j], linestyle = :solid, color = RGBAf(0.2, 0.5, 0.7, 0.5)) : _MK.Box(gl_matrix[i, j], linestyle = :dot, color = RGBAf((i*j)*0.01, 0.53, 0.93, 0.095))
                        zij_polar = polarize(zij)

                        Label(gl_matrix[i, j], "$(zij) \n $(zij_polar) \n ($(lc_imp[i, j]))", fontsize=14, tellwidth=false, tellheight=false, justification=:center, color=:black)
                    end
        end
    end


    # Visualize from bus
    visualize_bus_terminals(
        pf_sol["bus"][string(f_bus_math_id)]["voltage"], 
        string(f_bus_math_id), 
        branch_c_fr, 
        gl_fbus, 
        keep_pu, # keep_pu 
        Ibase_A, 
        vbase_V, 
        rounding_digits,
        string(math["bus"][string(f_bus_math_id)]["name"]) # pass the engineering bus id
    )

    # Visualize to bus
    visualize_bus_terminals(
        pf_sol["bus"][string(t_bus_math_id)]["voltage"], 
        string(t_bus_math_id), 
        branch_c_to, 
        gl_tbus,
        keep_pu, # keep_pu
        Ibase_A, 
        vbase_V, 
        rounding_digits,
        string(math["bus"][string(t_bus_math_id)]["name"]) # pass the engineering bus id

    )
    #_MK.Box(gl_fbus[1:end, 1], color = :transparent, strokecolor = :red, tellheight=true) # Left panel box
    #_MK.Box(gl_tbus[1:end, 1], color = :transparent, strokecolor = :red, tellheight=true) # Right panel box
    
    resize_to_layout!(f)
    return f


end



function visualize_bus_terminals(
    bus_Vs::Dict, 
    bus_id::String, 
    branch_currents::Dict, 
    gl_bus::GridLayout, 
    keep_pu::Bool=false, 
    Ibase_A::Real=1.0, 
    vbase_V::Real=1.0, 
    rounding_digits::Int=4,
    eng_id::String = "N/A" # Engineering bus id
)
    # Add title for the bus
    title_gl = GridLayout(gl_bus[0,1])
    Label(title_gl[1,1], "Bus $bus_id [$(eng_id)]", fontsize=18, tellwidth=true, justification=:center, color=:black)
    
    # Process each terminal
    for (t, voltage) in sort(bus_Vs)
        # Process current
        current = keep_pu ? branch_currents[t] : branch_currents[t] * Ibase_A
        current = round(current, digits=rounding_digits)
        current_polar = Pliers.polarize(current)
        
        # Process voltage
        voltage = keep_pu ? voltage : voltage * vbase_V
        voltage = round(voltage, digits=rounding_digits)
        voltage_polar = Pliers.polarize(voltage, rounding_digits=rounding_digits)
        
        # Get phase letter
        pl = Pliers._phase_letter(parse(Int, t))
        
        # Set up terminal grid layout
        terminal_gl = GridLayout(gl_bus[parse(Int, t), 1])
        
        # Terminal label
        Label(terminal_gl[1:4, 0], L"\textbf{%$(pl)}", fontsize=16, tellwidth=true, tellheight=true, justification=:center, color=:navy)
        
        # Voltages per terminal
        Label(terminal_gl[1:2, 1], L"U_%$(pl)", fontsize=16, tellwidth=true, tellheight=true, justification=:center, color=:royalblue)
        Label(terminal_gl[1, 2], "$(voltage)", fontsize=15, tellwidth=true, tellheight=true, justification=:center, color=:black)
        Label(terminal_gl[2, 2], "$(voltage_polar)", fontsize=15, tellwidth=true, tellheight=true, justification=:center, color=:black)
        
        # Currents per terminal
        Label(terminal_gl[3:4, 1], L"I_%$(pl)", fontsize=16, tellwidth=true, tellheight=true, justification=:center, color=:darkslateblue)
        Label(terminal_gl[3, 2], "$(current)", fontsize=15, tellwidth=true, tellheight=true, justification=:center, color=:black)
        Label(terminal_gl[4, 2], "$(current_polar)", fontsize=15, tellwidth=true, tellheight=true, justification=:center, color=:black)

        _MK.Box(terminal_gl[1:end, 0:2], color = :transparent, strokecolor = :black, tellheight=true, linestyle = :solid) 
        _MK.Box(terminal_gl[1:2, 1:2], color = :transparent, strokecolor = :royalblue, tellheight=true, linestyle = :dot) 
        _MK.Box(terminal_gl[3:4, 1:2], color = :transparent, strokecolor = :darkslateblue, tellheight=true, linestyle = :dot) 

    end
end



#=
 __     __         _______                       ______   __  __           
/  |   /  |       /       \                     /      \ /  |/  |          
$$ |   $$ |       $$$$$$$  | ______    ______  /$$$$$$  |$$/ $$ |  ______  
$$ |   $$ |______ $$ |__$$ |/      \  /      \ $$ |_ $$/ /  |$$ | /      \ 
$$  \ /$$//      |$$    $$//$$$$$$  |/$$$$$$  |$$   |    $$ |$$ |/$$$$$$  |
 $$  /$$/ $$$$$$/ $$$$$$$/ $$ |  $$/ $$ |  $$ |$$$$/     $$ |$$ |$$    $$ |
  $$ $$/          $$ |     $$ |      $$ \__$$ |$$ |      $$ |$$ |$$$$$$$$/ 
   $$$/           $$ |     $$ |      $$    $$/ $$ |      $$ |$$ |$$       |
    $/            $$/      $$/        $$$$$$/  $$/       $$/ $$/  $$$$$$$/ 
=#


## FIRST: GET LENGTHS
# I will try to create an aglorithm that starts from each node of the network, and then finds the line that connects to it and then goes to the parent node of that bus and so on unitl it reaches the root bus that has no lines going to it

# first a function that find the parent line of a bus
function find_parent_line(eng, bus_id)
    for (l, line) in eng["line"]
        if haskey(line, "t_bus") && line["t_bus"] == bus_id
            return l
        end
    end
end


# then a function that finds the parent bus of a line
function find_parent_bus(eng, line_id)
    return eng["line"][line_id]["f_bus"]
end


# define the root bus as the bus that doesn't have a line that takes it as a ["t_bus"]
function find_root_bus(eng)
    for (b, _) in eng["bus"]
        if isnothing(find_parent_line(eng, b))
            return b
        end
    end
end


# Now loop through all the buses and find the path to the root bus

function find_path_to_root(eng)
    root_bus = find_root_bus(eng)
    paths = Dict()
    for (b,_) in eng["bus"]
        papa = b
        paths[b] = []
        while papa != root_bus
            line_upstream = find_parent_line(eng, papa)
            papa = find_parent_bus(eng, line_upstream )
            push!(paths[b], line_upstream)
        end
    end
    return paths
end



# now I want a lengths dictionary for each bus that sums the length of the lines in the path

function find_path_lengths(eng)
    paths = find_path_to_root(eng)
    lengths = Dict()

    for (b, path) in paths
        lengths[b] = 0
        for l in path
            lengths[b] += eng["line"][l]["length"]
        end
    end
    return lengths
end


## SECOND: the plot of Voltages of a certain phase using bus["V"] value with length in X-axis and Voltage Magnitude in Y-axis



function distance_voltage_array(eng_res, Distances, phase)
    bus_labels = []
    x_distances_m = []
    y_voltages_V = []
    for (b, bus) in eng_res["bus"]
        Vbase = eng_res["bases"]["vbase_V"]
        display(b)
        if haskey(bus["V"], phase)
            push!(bus_labels, b)
            push!(x_distances_m, Distances[b])
            push!(y_voltages_V, abs(bus["V"][phase])*Vbase)
            println("Bus: $b, X: $(Distances[b]) Y= $(abs.(bus["V"][phase])*Vbase) ")
        end
    end
    return bus_labels, x_distances_m, y_voltages_V
end

function plot_voltage_profile(eng_res; size = (800, 1000), phase = nothing)
    @info "Calculating path lengths from each bus to the root bus"
    Distances = find_path_lengths(eng_res)
    @info "Phew done! d=====(￣▽￣*)b that took a while"

    f = Figure(size = size)
    colors = [:red, :green, :blue, :black]

    if isnothing(phase)
        axs = []
        for (i, color) in enumerate(colors)
            ax = Axis(f[i,1])
            bus_labels, x_distances_m, y_voltages_V = distance_voltage_array(eng_res, Distances, string(i))
            if isempty(bus_labels)
                display("No bus found for phase $( Dict("1" => "a", "2" => "b", "3" => "c", "4" => "n")[string(i)])") 
                return Distances
            end
            scatter!(ax, x_distances_m, y_voltages_V, color = color)
            ax.title = "Phase $( Dict("1" => "a", "2" => "b", "3" => "c", "4" => "n")[string(i)]) of the network"
            ax.xlabel = "Distance (m)"
            ax.ylabel = "Voltage (V)"    
            push!(axs, ax)
        end
        linkxaxes!(axs...)
        pop!(axs)
        display(axs)
        linkyaxes!(axs...)
    else
        ax = Axis(f[1,1])
        bus_labels, x_distances_m, y_voltages_V = distance_voltage_array(eng_res, Distances, string(phase))
        if isempty(bus_labels)
            display("No bus found for phase $( Dict("1" => "a", "2" => "b", "3" => "c", "4" => "n")[string(phase)])") 
            return Distances
        end
        scatter!(ax, x_distances_m, y_voltages_V, color = colors[phase])
        ax.title = "Phase $( Dict("1" => "a", "2" => "b", "3" => "c", "4" => "n")[string(phase)]) of the network"
        ax.xlabel = "Distance (m)"
        ax.ylabel = "Voltage (V)"
    end

    return f
end


# A phasor plotter

"""
    bus_phasor!(ax::PolarAxis, eng::Dict{String, Any}, bus_id::String;
                linestyle=:solid, colors=[:darkred, :darkgreen, :darkblue, :black])

Plots voltage phasors for a given bus onto an existing `PolarAxis`.

# Arguments
- `ax::PolarAxis`: The Makie PolarAxis to plot on.
- `eng::Dict{String, Any}`: The engineering data dictionary containing results.
- `bus_id::String`: The ID of the bus to plot phasors for.
- `linestyle`: The line style for the phasors (default: `:solid`).
- `colors`: The colors to use for the different phases (default: `[:darkred, :darkgreen, :darkblue, :black]`).
"""
function bus_phasor!(ax::PolarAxis, eng::Dict{String, Any}, bus_id::Integer;
                     linestyle=:solid, colors=[:darkred, :darkgreen, :darkblue, :black], keep_pu::Bool=false)

    # Extract base voltage - assuming results are not in per unit based on original logic
    is_perunit, vbase_V, sbase_VA, Zbase_Ω, Ibase_A, vbase_ll, Ibase_A_ll, Ibase_A_ϕ = calc_bases_from_dict(eng)
    bus_id = string(bus_id)
    bus = eng["bus"][bus_id]

    for (i, color) in enumerate(colors)
        phase_key = string(i)
        
        if haskey(bus["voltage"], phase_key)
            V_complex = bus["voltage"][phase_key]
            Vm = keep_pu ? abs(V_complex) : abs(V_complex) * vbase_V
            θ = angle(V_complex)
            lines!(ax, [0, θ], [0, Vm], color=color, linewidth=2, linestyle=linestyle)
        end
    end
    return ax
end


"""
    bus_phasor(eng::Dict{String, Any}, bus_id::String;
               makie_backend=WGLMakie, fig_size=(800, 800))

Creates a new figure and plots voltage phasors for a given bus.

# Arguments
- `eng::Dict{String, Any}`: The engineering data dictionary containing results.
- `bus_id::String`: The ID of the bus to plot phasors for.
- `makie_backend`: The Makie backend to activate (default: `WGLMakie`).
- `fig_size`: The size of the figure (default: `(800, 800)`).

# Returns
- `Figure`: The Makie figure containing the phasor plot.
"""
function bus_phasor(eng::Dict{String, Any}, bus_id::Integer;
                    makie_backend=WGLMakie,
                    figure::Figure = nothing,
                    location::Tuple{Int, Int} = (1, 1),
                    fig_size=(800, 800),
                    keep_pu::Bool=false,
                   )
    if isnothing(figure) 
        makie_backend.activate!()
        f = Figure(size=fig_size)
    else
        f = figure
    end

    degree_ticks = 0:30:330
    radian_ticks = deg2rad.(degree_ticks)
    tick_labels = ["$(d)°" for d in degree_ticks]

    ax = PolarAxis(f[location[1], location[2]],
                   title="Bus $bus_id Phasors",
                   thetaticks=(radian_ticks, tick_labels),
                  )

    bus_phasor!(ax, eng, bus_id, keep_pu=keep_pu)

    return f, ax
end

function Vphasor(data::Dict{String, Any}, bus_id::Integer;
         makie_backend=WGLMakie, fig_size=(800, 800),colors=[:darkred, :darkgreen, :darkblue, :black],keep_pu = true, kwargs...) 
            
         makie_backend.activate!()

         f = Figure(size=fig_size)
            degree_ticks = 0:30:330
            radian_ticks = deg2rad.(degree_ticks)
            tick_labels = ["$(d)°" for d in degree_ticks]

            ax = PolarAxis(f[1,1],
                           title = "Bus $bus_id Voltage Phasors",
                           thetaticks = (radian_ticks, tick_labels),
                           kwargs...
                          )
        is_perunit, vbase_V, sbase_VA, Zbase_Ω, Ibase_A, vbase_ll, Ibase_A_ll, Ibase_A_ϕ = calc_bases_from_dict(data)
        bus_id = string(bus_id)
        bus = data["bus"][bus_id]

        for (i, color) in enumerate(colors)
            phase_key = string(i)
            
            if haskey(bus["voltage"], phase_key)
                V_complex = bus["voltage"][phase_key]
                Vm = keep_pu ? abs(V_complex) : abs(V_complex) * vbase_V
                θ = angle(V_complex)
                lines!(ax, [0, θ], [0, Vm], color=color, linewidth=2, linestyle=:solid)
            end
        end

        return f
end


"""
    vphasor(scene, value; kwargs...)
    vphasor!(scene, value; kwargs...)

Plot a complex number as a phasor (arrow from origin) in a 2D plane.

# Arguments
- `scene`: The Makie scene to plot into
- `value`: A complex number to be plotted as a phasor

# Attributes
- `arrowcolor`: Color of the phasor arrow (default: current theme's linecolor)
- `arrowsize`: Size of the arrowhead in pixels (default: 15)
- `linewidth`: Width of the phasor line (default: current theme's linewidth)
- `label`: Label for the phasor in the legend (default: "")
- `linestyle`: Style of the phasor line (default: nothing)
- `inspectable`: Whether the phasor is inspectable with mouse hover (default: true)

# Examples
v1 = 1 + 2im
v2 = -2 - 1im
v3 = 0.5 - 1.5im

vphasor(scene, v1, arrowcolor = :red, label = "V1")
vphasor!(scene, v2, arrowcolor = :green, label = "V2")
vphasor!(scene, v3, arrowcolor = :blue, label = "V3")
"""
# VPhasor recipe: plots a complex number as an arrow from origin.
MakieCore.@recipe(VPhasor, value) do scene
    Theme(
        arrowcolor = Makie.theme(scene, :linecolor), # Default to current linecolor from the theme
        arrowsize = 15,
        linewidth = Makie.theme(scene, :linewidth), # Default to current linewidth from the theme
        label = "",
        linestyle = nothing, # Allow linestyle to be passed to arrows!
        inspectable = true   # Default for inspectable attribute
    )
end

function Makie.plot!(plot::VPhasor{<:Tuple{C}}) where {C <: Complex}
    # value_obs is an Observable{Complex} containing the input complex number
    value_obs = plot.value

    # Create an Observable for the arrow points
    # The arrow will go from the origin (0,0) to (real(value), imag(value))
    arrow_points = lift(value_obs; ignore_equal_values=true) do val
        [Point2f(0, 0), Point2f(real(val), imag(val))]
    end

    # Use Makie.arrows! to draw the phasor
    # plot.attributes gives access to the theme and any attributes passed by the user
    Makie.arrows!(plot, arrow_points;
        color = plot.arrowcolor,
        arrowsize = plot.arrowsize,
        linewidth = plot.linewidth,
        linestyle = plot.linestyle,
        label = plot.label,
        inspectable = plot.inspectable
    )

    return plot # Return the plot object
end





function line_current_phasor(eng::Dict{String, Any}, line_id::String; keep_pu::Bool= false, makie_backend=WGLMakie,
    )
    Ibase_A = eng["bases"]["Ibase_A"]

    makie_backend.activate!()

    f = Figure(size=(800, 1200))

    colors = [:red, :green, :blue, :black]

    degree_ticks = 0:30:330
    radian_ticks = deg2rad.(degree_ticks)
    tick_labels = ["$(d)°" for d in degree_ticks]

    ax = PolarAxis(f[1,1],
                   title = "Line $line_id Current Phasors (From Bus)",
                   thetaticks = (radian_ticks, tick_labels)
                  )

    line = eng["line"][line_id]
    for (i, color) in enumerate(colors)
        phase_key = string(i)
        if haskey(line["I_f"], phase_key)
            I_complex = line["I_f"][phase_key]
            Im = abs(I_complex) * Ibase_A
            θ = angle(I_complex)
            lines!(ax, [0,θ], [0,Im], color = color, linewidth = 2, linestyle = :solid)
        end
    end

    return f 
end