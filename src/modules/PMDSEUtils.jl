"""
    PMDSEUtils

Internal sub-module providing utility functions for PowerModelsDistributionStateEstimation (PMDSE) workflows.

This module re-exports functions from the main Pliers module for:
- State estimation result visualization
- Measurement residual analysis
- Measurement data processing and writing

See the main Pliers module for function documentation.
"""
module PMDSEUtils

using ..Pliers
using ..Pliers: warning_text, header, extra_keys
using ..Pliers.PMDUtils: dictify_solution!, add_vmn_p_q
using ..Pliers.PMDGraph: create_network_graph, smart_layout,
    network_graph_plot, _decorate_nodes!, _decorate_edges!

# plotting packages
using Makie
using CairoMakie
using WGLMakie
using Graphs, MetaGraphs

# pretty terminal packages
using Crayons
using Crayons.Box
using PrettyTables

using DataFrames
using CSV

# data analysis packages
using Statistics
using LinearAlgebra


_N_IDX = 4 # representing the neutral index in terminals 


function get_m_n_dof(data::Dict)

    ref_bus = [bus for (_,bus) in data["bus"] if bus["bus_type"] == 3]
    
    @assert length(ref_bus) == 1 "There is more than one reference bus, double-check model"
    
    load_buses = ["$(load["load_bus"])" for (_, load) in data["load"]] # buses with demand (incl. negative demand, i.e., generation passed as negative load)
    gen_slack_buses = ["$(gen["gen_bus"])" for (_, gen) in data["gen"]] # buses with generators, including the slackbus
    non_zero_inj_buses = unique(vcat(load_buses, gen_slack_buses))
    
    @assert !isempty(non_zero_inj_buses) "This network has no active connected component, no point doing state estimation"
    
    n = sum([length(setdiff(bus["terminals"],_N_IDX)) for (b, bus) in data["bus"] if b ∈ non_zero_inj_buses ])*2-length(setdiff(ref_bus[1]["terminals"],_N_IDX)) # system variables: two variables per bus (e.g., angle and magnitude) per number of phases on the bus minus the angle variable(s) of a ref bus, which are fixed
    #n = sum([length(bus["terminals"]) for (b, bus) in data["bus"] if b ∈ non_zero_inj_buses ])*2-length(ref_bus[1]["terminals"]) # system variables: two variables per bus (e.g., angle and magnitude) per number of phases on the bus minus the angle variable(s) of a ref bus, which are fixed
    #@show n
    m = sum([length(meas["dst"]) for (_, meas) in data["meas"]])
    #@show m
    m-n > 0 ? nothing : warning_text("system underdetermined or just barely determined")
    
    return m, n, m-n

end


#TODO: fix it to not only show 3 columns but extend the res array,, also note math has maybe all things u need?!
function viz_residuals(SE_en, math_en; show_legend = false, mape = nothing, MAPE_N = nothing, APEs_df = nothing)#pm_form = PowerModelsDistribution.IVRENPowerModel)

    se_sol_en = SE_en["solution"]
    solve_time = SE_en["solve_time"]
    objective = SE_en["objective"]
    termination_status = SE_en["termination_status"]
    primal_status = SE_en["primal_status"]
    

    printstyled(" %%%%%%%%%%%%%%%%%% STATS %%%%%%%%%%%%%%%%%% \n", color=:blue, underline = true)
    printstyled(" Solve time : $(solve_time) \n", color=:green)
    isapprox(objective, 0, atol = 0.05) ? printstyled(" objective : $(objective) \n", color=:green) : printstyled(" objective : $(objective) \n", color=:red)
    string(termination_status) == "LOCALLY_SOLVED" ? printstyled(" termination : $(termination_status) \n", color=:green) : printstyled(" termination : $(termination_status) \n", color=:red)
    string(primal_status) == "FEASIBLE_POINT" ? printstyled(" primal : $(primal_status) \n", color=:green) : printstyled(" primal : $(primal_status) \n", color=:red)
    !isnothing(mape) ? isapprox(mape, 0, atol=0.05)  ? printstyled(" MAPE : $(mape) \n", color=:green) : printstyled(" MAPE : $(mape) \n", color=:red) : nothing
    !isnothing(MAPE_N) ? isapprox(MAPE_N, 0, atol=0.05)  ? printstyled(" MAPE_N : $(MAPE_N) \n", color=:green) : printstyled(" MAPE_N : $(MAPE_N) \n", color=:red) : nothing
    
    m,n, dof = get_m_n_dof(math_en)
    printstyled(" m : $(m) \n", color=:green)
    printstyled(" n : $(n) \n", color=:green)
    dof > 0 ? printstyled(" Degrees of freedom : $(dof) \n", color=:green) : printstyled(" Degrees of freedom : $(dof) \n", color=:red)
    
    # @suppress display(instantiate_mc_model(
    #     math_en,
    #     pm_form,
    #     PowerModelsDistributionStateEstimation.build_mc_se;
    # ).model)

    if !isapprox(mape, 0, atol=0.05)
        if !isnothing(APEs_df)
            println("APEs_df:")
            display(APEs_df)
        end
    end

    # Merge residuals with meas
    for (m, meas) in se_sol_en["meas"]
        math_en["meas"][m]["res"] = meas["res"]
    end

    # Extracting values into a DataFrame for visualization
    rows = String[]
    res1 = Float64[]
    res2 = Float64[]
    res3 = Float64[]
    res4 = Float64[]

    for (key, value) in math_en["meas"]
        cmp = value["cmp"]
        cmp_id = value["cmp_id"]
        var = value["var"]
        name = math_en[string(cmp)][string(cmp_id)]["name"]
        row_name = "$(cmp)_$(cmp_id)_$(var)"
        push!(rows, row_name)
        res_values = value["res"]
        push!(res1, get(res_values, 1, NaN))
        push!(res2, get(res_values, 2, NaN))
        push!(res3, get(res_values, 3, NaN))
        push!(res4, get(res_values, 4, NaN))
    end

    df_meas_res = DataFrame(Meas = rows, Res1 = res1, Res2 = res2, Res3 = res3, Res4 = res4)
    sort!(df_meas_res, :Meas)

    # Handle empty collections for minimum and maximum calculations
    non_nan_res1 = filter(!isnan, df_meas_res.Res1)
    non_nan_res2 = filter(!isnan, df_meas_res.Res2)
    non_nan_res3 = filter(!isnan, df_meas_res.Res3)
    non_nan_res4 = filter(!isnan, df_meas_res.Res4)

    min_val = minimum(vcat(non_nan_res1, non_nan_res2, non_nan_res3, non_nan_res4))
    max_val = maximum(vcat(non_nan_res1, non_nan_res2, non_nan_res3, non_nan_res4))

    # Define color intervals
    num_intervals = 6
    lowerBounds = [[min_val + i * (0.01 - min_val) / 3 for i in 0:3]..., 1e-1, 1]
    uppoerBounds = [[0.00333332833981672 + i * (0.1 - 0.00333332833981672) / 3 for i in 0:3]..., 1, max_val + 1e-6]

    colors = [crayon"green bold", crayon"cyan bold", crayon"blue bold", crayon"magenta bold", crayon"yellow bold", crayon"red bold"]

    # Create highlighters for each interval
    highlighters = TextHighlighter[]
    for i in 1:num_intervals
        lower_bound = lowerBounds[i]
        upper_bound = uppoerBounds[i]
        push!(highlighters, TextHighlighter(
            (data, row, col) -> col in 2:5 && data[row, col] >= lower_bound && data[row, col] < upper_bound,
            colors[i]
        ))
    end
    hideNan_hl = TextHighlighter(       
        (data, i, j) -> (j ∈ collect(2:5) && isnan(data[i, j])),
        crayon"dark_gray conceal"
    )
    push!(highlighters, hideNan_hl)

    # Define headers
    header = ["Meas", "Res1", "Res2", "Res3", "Res4"]
    try
        res_heatmap = pretty_table(df_meas_res; column_labels=header, highlighters=highlighters)
    catch
        display(df_meas_res)
    end
    df_legend = DataFrame(
        Lower_Bound = lowerBounds,
        Upper_Bound = uppoerBounds,
    )

    legend_highlighters = TextHighlighter[]
    for i in 1:num_intervals
        push!(legend_highlighters, TextHighlighter(
            (data, row, col) -> (row == i),
            colors[i]
        ))
    end
    if show_legend
        println("Legend:") 
        legend = pretty_table(df_legend, header=["Lower Bound", "Upper Bound"], header_crayon=crayon"fg:yellow", highlighters=legend_highlighters)
        println(legend)
    end


    return solve_time, objective, termination_status, primal_status, mape, m, n, dof

end

function viz_residuals(SE_en, math_en, PF_en; show_legend=false, kwargs...)

    MAPE, APEs_df, _ = _calculate_MAPE(SE_en, PF_en, math_en)
    MAPE_N, _ , _ = _calculate_MAPE_toNeutral(SE_en, PF_en, math_en)
    viz_residuals(SE_en, math_en; show_legend = show_legend, mape = MAPE, APEs_df = APEs_df, MAPE_N=MAPE_N, kwargs...)

end



function df_meas_res(SE_en, math_en)

    se_sol_en = SE_en["solution"]
        
    printstyled(" %%%%%%%%%%%%%%%%%% STATS %%%%%%%%%%%%%%%%%% \n", color=:blue, underline = true)
    printstyled(" Solve time : $(SE_en["solve_time"]) \n", color=:green)
    printstyled(" objective : $(SE_en["objective"]) \n", color=:green)
    printstyled(" termination : $(SE_en["termination_status"]) \n", color=:green)
    printstyled(" primal : $(SE_en["primal_status"]) \n", color=:green)

    # Merge residuals with meas
    for (m, meas) in se_sol_en["meas"]
        math_en["meas"][m]["res"] = meas["res"]
    end


    # Extracting values into a DataFrame for visualization
    rows = String[]
    res1 = Float64[]
    res2 = Float64[]
    res3 = Float64[]

    for (key, value) in math_en["meas"]
        cmp = value["cmp"]
        cmp_id = value["cmp_id"]
        var = value["var"]
        name = math_en[string(cmp)][string(cmp_id)]["name"]
        row_name = "$(cmp)_$(cmp_id)_$(var)"
        push!(rows, row_name)
        res_values = value["res"]
        push!(res1, get(res_values, 1, NaN))
        push!(res2, get(res_values, 2, NaN))
        push!(res3, get(res_values, 3, NaN))
    end

    df_meas_res = DataFrame(Meas = rows, Res1 = res1, Res2 = res2, Res3 = res3)
    sort!(df_meas_res, :Meas)

    return df_meas_res
end


# Residuals tabls for SE solutions 


function math_meas_table(math::Dict{String, Any}; se_sol= nothing)
    

    meass = haskey(math, "meas") ?  math["meas"] : error("No measurements found in the math data models")

    # Merge residuals with meas
    if !isnothing(se_sol)
        if haskey(se_sol, "meas")
            for (m, meas) in math["meas"]
                meas["res"] = se_sol["meas"][m]
            end
        else
            warning_text("No measurements residuals `se_sol[\"meas\"]` found in the SE solution")
            for (_, meas) in math["meas"]
                meas["res"] = nothing
            end
        end 
    else
        if haskey(math, "res") 
            
        else
            warning_text("No measurements residuals `math[\"res\"]` found in the math data models")
            for (_, meas) in math["meas"]
                meas["res"] = nothing 
            end
        end 
    end
    
    # iterate through all math["meas"] and write a table 

    #=
    math_1ph["meas"]["10"]
    Dict{String, Any} with 6 entries:
    "var"    => :cid
    "crit"   => "rwlav"
    "cmp"    => :load
    "res"    => [0.432542, -7.4941e-9]
    "dst"    => Normal{Float64}[Distributions.Normal{Float64}(μ=7.39022e-7, σ=0.00333333), Distributions.Normal{Float64}(μ=-7.39022e-7, σ=0.00333333)]
    "cmp_id" => 29
    =#


    

    meas_df = DataFrame(Components = String[], Component_ID = Int[], Variable = String[], Criterion = String[], Residuals = Union{Vector{Float64}, Nothing}[])
    for (_, value) in meass
            cmp = value["cmp"]
            cmp_id = value["cmp_id"]
            var = value["var"]
            crit = value["crit"]
            name = math[string(cmp)][string(cmp_id)]["name"]
            push!(meas_df, (string(cmp), cmp_id, string(var), string(crit), value["res"]))            
        end
        header("Measurements Table: contains ($(nrow(meas_df)) measurements)")
        
        headers = (
           ["Components"  "Component_ID"  "Variable"  "Criterion "  "Residuals" ],
           [ "cmp"  "cmp_id"  "var"  "crit "  "res"]
        )
        pretty_table(meas_df, header= headers)
        extra_keys(meass, ["cmp", "cmp_id", "var", "crit", "res"])
    
end

#TODO: refactor to be more efficient
function math_meas_table(math::Dict{String, Any}, condition::Function; se_sol= nothing)
    meass = haskey(math, "meas") ?  math["meas"] : error("No measurements found in the math data models")
    # Merge residuals with meas
    if !isnothing(se_sol)
        if haskey(se_sol, "meas")
            for (m, meas) in math["meas"]
                meas["res"] = se_sol["meas"][m]
            end
        else
            warning_text("No measurements residuals `se_sol[\"meas\"]` found in the SE solution")
            for (_, meas) in math["meas"]
                meas["res"] = nothing
            end
        end 
    else
        if haskey(math, "res") 
            
        else
            warning_text("No measurements residuals `math[\"res\"]` found in the math data models")
            for (_, meas) in math["meas"]
                meas["res"] = nothing 
            end
        end 
    end



    meas_df = DataFrame(Components = String[], Component_ID = Int[], Variable = String[], Criterion = String[], Residuals = Union{Vector{Float64}, Nothing}[])
    for (_, value) in meass
            cmp = value["cmp"]
            cmp_id = value["cmp_id"]
            var = value["var"]
            crit = value["crit"]
            name = math[string(cmp)][string(cmp_id)]["name"]
            if condition(value)
                push!(meas_df, (string(cmp), cmp_id, string(var), string(crit), value["res"]))            
            end            
        end
    header("Measurements Table: contains ($(nrow(meas_df)) measurements)")
    headers = (
        ["Components"  "Component_ID"  "Variable"  "Criterion "  "Residuals" ],
        [ "cmp"  "cmp_id"  "var"  "crit "  "res"]
     )
     pretty_table(meas_df, header= headers)
    extra_keys(meass, ["cmp", "cmp_id", "var", "crit", "res"])#
        
end


function write_sm_measurements(PF_RES, math, measurements_file, measurement_model, write_measurements_fn::Function; σ=0.05)
    dictify_solution!(PF_RES["solution"], math)
    math_meas_en = add_vmn_p_q(math, PF_RES["solution"])
    write_measurements_fn(measurement_model, math_meas_en, PF_RES, measurements_file, σ=σ) 
end


function add_pd_qd_vmn!(SE_RES::Dict{String, Any}, math::Dict{String, Any})
    for (l, lo) in SE_RES["solution"]["load"]
        load_bus = math["load"][l]["load_bus"]
        
        lo["power"] = Dict{}()
        for (t, I_l_cmplx) in lo["current"]
            V_l_cmplx = SE_RES["solution"]["bus"][string(load_bus)]["voltage"][t]
            S_l_cmplx = V_l_cmplx * conj(I_l_cmplx)
            lo["power"][t] = Dict("pd" => real(S_l_cmplx), "qd" => imag(S_l_cmplx), "S" => S_l_cmplx)
        end
        sorted_keys = sort(collect(keys(lo["power"])), by = x -> parse(Int, x))
        lo["pd"] = [lo["power"][key]["pd"] for key in sorted_keys]
        lo["qd"] = [lo["power"][key]["qd"] for key in sorted_keys]
    end
    
    for (b, bus) in SE_RES["solution"]["bus"]
        bus["vmn"] = Dict{}()
        for (t, v) in bus["voltage"]
            if t in ["1", "2", "3"]
                bus["vmn"][t] = haskey( bus["voltage"], "4") ? abs(v - bus["voltage"]["4"]) : abs(v)
            end
        end
        sorted_keys = sort(collect(keys(bus["vmn"])), by = x -> parse(Int, x))
        bus["vmn"] = [bus["vmn"][key] for key in sorted_keys]
    end
end  



"""
    _calculate_MAPE(SE_RES, PF_RES, math)

Calculate the Mean Absolute Percentage Error (MAPE) between state estimation (SE) results 
and power flow (PF) results.

# Arguments
- `SE_RES::Dict`: A dictionary containing the state estimation results.
- `PF_RES::Dict`: A dictionary containing the power flow results.
- `math`: A mathematical utility or context used for processing the solutions.

# Returns
- `mean_APE::Float64`: The mean absolute percentage error across all buses and terminals.
- `APEs::Vector{Float64}`: A vector containing the absolute percentage errors for each bus and terminal.

# Details
The function compares the voltage solutions from the state estimation (`SE_RES`) and 
power flow (`PF_RES`) results. It calculates the absolute percentage error (APE) for 
each bus and terminal, filters out any `NaN` values, and computes the mean APE.

# Notes
- The function assumes that the solutions in `SE_RES` and `PF_RES` are structured as 
  dictionaries with nested keys for buses and their respective voltages.
- The `dictify_solution!` function is used to process the solutions before comparison.
"""
function _calculate_MAPE(SE_RES, PF_RES, math)

    se_sol  = deepcopy(SE_RES["solution"])
    pf_sol = deepcopy(PF_RES["solution"])

    dictify_solution!(pf_sol, math)
    dictify_solution!(se_sol, math)

    APEs = [] 
    Errors = [] 

    Errors_df = DataFrame(Bus = String[], Terminal = String[], Error = ComplexF64[])
    APEs_df = DataFrame(Bus = String[], Terminal = String[], APE = Float64[])
    for (b,bus) in se_sol["bus"]

        for (term,Vse) in bus["voltage"] 
            if term == "4"
                continue
            end
            # println("Bus    :",b)
            # println("Term   :",term)
            # println("Vse   :",Vse)
            # println("Vpf     :", pf_sol["bus"][b]["voltage"][term])


            #haskey(bus["voltage"], "4") ? nothing : bus["voltage"]["4"] = 0.0 + 0.0im # add a dummy voltage for the 4th terminal if it doesn't exist
            #haskey(pf_sol["bus"][b]["voltage"], "4") ? nothing : bus["voltage"]["4"] = 0.0 + 0.0im # add a dummy voltage for the 4th terminal if it doesn't exist
            Vpf = pf_sol["bus"][b]["voltage"][term] # subtract the dummy voltage for the 4th terminal if it doesn't exist
            #Vpf = haskey(pf_sol["bus"][b]["voltage"], "4")  ?  pf_sol["bus"][b]["voltage"][term] -  pf_sol["bus"][b]["voltage"]["4"]   : pf_sol["bus"][b]["voltage"][term] # subtract the dummy voltage for the 4th terminal if it doesn't exist
            #Vse = haskey(bus["voltage"] , "4") ? Vse - bus["voltage"]["4"] : Vse
            # APE = abs.(Vpf) == 0 ? 0 : abs.( Vse - Vpf ) ./ abs.(Vpf) * 100
            #APE = abs.(abs.( Vse) .- abs.(Vpf)) ./ abs.(Vpf)* 100
            #APE = abs.(abs.( Vse) .- abs.(Vpf)) ./ abs.(Vpf)* 100
            APE = (abs.(Vse .- Vpf) ./ abs.(Vpf) ) * 100
            Error = Vse - Vpf
            push!(APEs,APE)
            push!(Errors,Error)


            push!(APEs_df, (b, term, APE))
            push!(Errors_df, (b, term, Error))

        end

    end
    
    APEs_nonan = filter(x -> !isnan(x), APEs)
    mean_APE = mean(APEs_nonan)
    return mean_APE, APEs_df, Errors_df

end

function _calculate_MAPE_toNeutral(SE_RES, PF_RES, math)

    se_sol  = deepcopy(SE_RES["solution"])
    pf_sol = deepcopy(PF_RES["solution"])

    dictify_solution!(pf_sol, math)
    dictify_solution!(se_sol, math)

    APEs = [] 
    Errors = [] 

    Errors_df = DataFrame(Bus = String[], Terminal = String[], Error = ComplexF64[])
    APEs_df = DataFrame(Bus = String[], Terminal = String[], APE = Float64[])
    for (b,bus) in se_sol["bus"]

        for (term,Vse) in bus["voltage"] 
            if term == "4"
                continue
            end
            # println("Bus    :",b)
            # println("Term   :",term)
            # println("Vse   :",Vse)
            # println("Vpf     :", pf_sol["bus"][b]["voltage"][term])


            #haskey(bus["voltage"], "4") ? nothing : bus["voltage"]["4"] = 0.0 + 0.0im # add a dummy voltage for the 4th terminal if it doesn't exist
            #haskey(pf_sol["bus"][b]["voltage"], "4") ? nothing : bus["voltage"]["4"] = 0.0 + 0.0im # add a dummy voltage for the 4th terminal if it doesn't exist

            #Vpf = pf_sol["bus"][b]["voltage"][term] 
            Vpf = haskey(pf_sol["bus"][b]["voltage"], "4")  ?  pf_sol["bus"][b]["voltage"][term] -  pf_sol["bus"][b]["voltage"]["4"]   : pf_sol["bus"][b]["voltage"][term] # subtract the dummy voltage for the 4th terminal if it doesn't exist
            Vse = haskey(bus["voltage"] , "4") ? Vse - bus["voltage"]["4"] : Vse
            # APE = abs.(Vpf) == 0 ? 0 : abs.( Vse - Vpf ) ./ abs.(Vpf) * 100
            #APE = abs.(abs.( Vse) .- abs.(Vpf)) ./ abs.(Vpf)* 100
            #APE = abs.(abs.( Vse) .- abs.(Vpf)) ./ abs.(Vpf)* 100
            APE = (abs.(Vse .- Vpf) ./ abs.(Vpf) ) * 100
            Error = Vse - Vpf
            push!(APEs,APE)
            push!(Errors,Error)


            push!(APEs_df, (b, term, APE))
            push!(Errors_df, (b, term, Error))

        end

    end
    
    APEs_nonan = filter(x -> !isnan(x), APEs)
    mean_APE = mean(APEs_nonan)
    return mean_APE, APEs_df, Errors_df

end

function _calculate_MAPE_P_Q(SE_RES, PF_RES, math; element = "load", quantity = "power")

    se_sol  = deepcopy(SE_RES["solution"])
    pf_sol = deepcopy(PF_RES["solution"])

    dictify_solution!(pf_sol, math)
    dictify_solution!(se_sol, math)

    P_APEs = [] 
    Q_APEs = []
    Errors = [] 

    Errors_df = DataFrame(Bus = String[], Terminal = String[], Error = ComplexF64[])
    P_APEs_df = DataFrame(Bus = String[], Terminal = String[], APE = Float64[])
    Q_APEs_df = DataFrame(Bus = String[], Terminal = String[], APE = Float64[])
    for (b,Component) in se_sol[element]

        for (term,Sse) in Component[quantity] 
            if term == "4"
                continue
            end

            Spf = pf_sol[element][b][quantity*"_bus"][term] # subtract the dummy voltage for the 4th terminal if it doesn't exist

            Ppf = real(Spf)
            Qpf = imag(Spf)

            Pse = real(Sse)
            Qse = imag(Sse)

            P_APE = (abs(abs(Pse) - abs(Ppf)) / abs(Ppf) ) * 100
            if P_APE > 100
                warning_text("High P_APE detected at Bus: $b, Terminal: $term, Pse: $Pse, Ppf: $Ppf, Sse: $Sse, Spf: $Spf")
            end
            Q_APE = (abs(abs(Qse) - abs(Qpf)) / abs(Qpf) ) * 100

            Error = Sse - Spf
            
            push!(P_APEs,P_APE)
            push!(Q_APEs,Q_APE)

            push!(Errors,Error)


            push!(P_APEs_df, (b, term, P_APE))
            push!(Q_APEs_df, (b, term, Q_APE))
            push!(Errors_df, (b, term, Error))

        end

    end

    P_APEs_nonan = filter(x -> !isnan(x), P_APEs)
    Q_APEs_nonan = filter(x -> !isnan(x), Q_APEs)
    
    mean_P_APE = mean(P_APEs_nonan)
    mean_Q_APE = mean(Q_APEs_nonan)


    return mean_P_APE, P_APEs_df, mean_Q_APE, Q_APEs_df, Errors_df

end


"""
    generate_meas_est_table(se_results, math)

Generates a DataFrame comparing measured values against estimated values from state estimation results.

# Arguments
- `se_results`: A dictionary containing the results of the state estimation, specifically looking for the `"solution"` key to access estimated variables.
- `math`: A dictionary containing the mathematical model and measurement data, specifically looking for the `"meas"` key which iterates over measurement definitions.

# Returns
- `DataFrame`: A DataFrame with the following columns:
    - `MeasID`: The unique identifier for the measurement.
    - `CompID`: The component identifier, concatenated as "component_type.component_id".
    - `Variable`: The specific variable being measured (e.g., voltage, power).
    - `MeasValue`: The mean value of the measurement distribution.
    - `EstValue`: The estimated value from the state estimation solution corresponding to the specific component and variable.

# Example
    df_results = generate_meas_est_table(se_results, math)

"""

function generate_meas_est_table(se_results, math)
    df_meas_est = DataFrame(MeasID=String[], CompID=String[], Variable=String[], MeasValue=Vector{Float64}[], EstValue=Vector{Float64}[])

    for (meas_id, meas) in sort(math["meas"])
        push!(df_meas_est, (MeasID=meas_id, CompID=string(meas["cmp"]) * "." * string(meas["cmp_id"]), Variable=string(meas["var"]), MeasValue=mean.(meas["dst"]), EstValue=se_results["solution"][string(meas["cmp"])][string(meas["cmp_id"])][string(meas["var"])]))
    end

    return df_meas_est
end

















# ═══════════════════════════════════════════════════════════════════════════════
# Residual Overlay Visualization on Network Graph
# ═══════════════════════════════════════════════════════════════════════════════

# Variable color and label mappings for residual bar charts
const _VAR_COLORS = Dict{Symbol, Symbol}(
    :vr  => :steelblue,
    :vi  => :coral,
    :pd  => :forestgreen,
    :qd  => :darkorange,
    :cid => :mediumpurple,
    :vmn => :gold,
    :pg  => :teal,
    :qg  => :salmon,
    :cr  => :royalblue,
    :ci  => :tomato,
)

const _VAR_LABELS = Dict{Symbol, String}(
    :vr  => "Vr",
    :vi  => "Vi",
    :pd  => "Pd",
    :qd  => "Qd",
    :cid => "Cid",
    :vmn => "Vmn",
    :pg  => "Pg",
    :qg  => "Qg",
    :cr  => "Cr",
    :ci  => "Ci",
)

const _VAR_ORDER = Dict{Symbol, Int}(
    :vr => 1, :vi => 2, :pd => 3, :qd => 4,
    :cid => 5, :vmn => 6, :pg => 7, :qg => 8, :cr => 9, :ci => 10,
)


"""
    _collect_bus_residuals(SE_en, math)

Aggregate measurement residuals by bus. Maps each measurement to its associated bus
based on the component type (`:bus`, `:load`, `:gen`).

Returns a `Dict{String, Vector}` where each key is a bus_id string and each value
is a vector of named tuples `(var, cmp, cmp_id, res)`.
"""
function _collect_bus_residuals(SE_en, math)
    # Merge residuals from solution into math
    for (m, meas) in SE_en["solution"]["meas"]
        math["meas"][m]["res"] = meas["res"]
    end

    bus_residuals = Dict{String, Vector{Any}}()

    for (_, meas) in math["meas"]
        cmp = meas["cmp"]
        cmp_id = meas["cmp_id"]
        var = meas["var"]
        res = meas["res"]

        bus_id = if cmp == :bus
            string(cmp_id)
        elseif cmp == :load
            string(math["load"][string(cmp_id)]["load_bus"])
        elseif cmp == :gen
            string(math["gen"][string(cmp_id)]["gen_bus"])
        else
            nothing
        end

        isnothing(bus_id) && continue

        if !haskey(bus_residuals, bus_id)
            bus_residuals[bus_id] = []
        end
        push!(bus_residuals[bus_id], (var=var, cmp=cmp, cmp_id=cmp_id, res=res))
    end

    return bus_residuals
end


"""
    _get_node_positions(network_graph, layout_fn)

Compute node positions from the graph and layout function.
Returns a `Dict{String, Tuple{Float64,Float64}}` mapping bus_id to `(x, y)`.
"""
function _get_node_positions(network_graph, layout_fn)
    positions = layout_fn(network_graph)
    bus_positions = Dict{String, Tuple{Float64, Float64}}()
    for v in vertices(network_graph)
        bus_id = string(get_prop(network_graph, v, :bus_id))
        bus_positions[bus_id] = positions[v]
    end
    return bus_positions
end


"""
    _aggregate_bus_residuals(bus_residuals, phase, aggregation)

Aggregate residual values per bus into a single scalar for bubble/heatmap visualizations.
Returns `Dict{String, Float64}`.
"""
function _aggregate_bus_residuals(bus_residuals, phase::Int, aggregation::Symbol)
    bus_agg = Dict{String, Float64}()
    for (bus_id, meas_list) in bus_residuals
        vals = Float64[]
        for m in meas_list
            val = get(m.res, phase, NaN)
            !isnan(val) && push!(vals, abs(val))
        end
        if !isempty(vals)
            bus_agg[bus_id] = if aggregation == :max
                maximum(vals)
            elseif aggregation == :mean
                mean(vals)
            elseif aggregation == :sum
                sum(vals)
            else
                maximum(vals)
            end
        end
    end
    return bus_agg
end


"""
    _base_network_plot(math, network_graph; kwargs...)

Create the base network topology plot shared by all residual visualization styles.
Returns `(Figure, Axis, GraphPlot)`.
"""
function _base_network_plot(math, network_graph;
    layout=smart_layout,
    show_node_labels::Bool=false,
    show_edge_labels::Bool=false,
    figure_size=(1000, 1200),
    makie_backend=WGLMakie,
    node_color_override=nothing,
    node_size_override=nothing,
    kwargs...)

    nlabels = show_node_labels ? [string(props(network_graph, i)[:bus_id]) for i in 1:nv(network_graph)] : nothing

    _decorate_nodes!(network_graph, math)
    node_color = isnothing(node_color_override) ?
        [props(network_graph, i)[:node_color] for i in 1:nv(network_graph)] : node_color_override
    node_marker = [props(network_graph, i)[:node_marker] for i in 1:nv(network_graph)]
    node_size = isnothing(node_size_override) ?
        [props(network_graph, i)[:marker_size] for i in 1:nv(network_graph)] : node_size_override

    _decorate_edges!(network_graph, math)
    edge_color = [get_prop(network_graph, e, :edge_color) for e in edges(network_graph)]
    arrow_show = true
    arrow_marker = [get_prop(network_graph, e, :arrow_marker) for e in edges(network_graph)]
    arrow_size = [get_prop(network_graph, e, :arrow_size) for e in edges(network_graph)]
    arrow_shift = [get_prop(network_graph, e, :arrow_shift) for e in edges(network_graph)]

    fig_ax_plot = network_graph_plot(
        network_graph;
        layout=layout,
        figure_size=figure_size,
        makie_backend=makie_backend,
        show_node_labels=show_node_labels,
        nlabels=nlabels,
        show_edge_labels=show_edge_labels,
        node_color=node_color,
        node_marker=node_marker,
        node_size=node_size,
        edge_color=edge_color,
        arrow_show=arrow_show,
        arrow_marker=arrow_marker,
        arrow_size=arrow_size,
        arrow_shift=arrow_shift,
        kwargs...
    )

    return fig_ax_plot
end


"""
    plot_residuals_bars(SE_en, math; kwargs...)

Plot measurement residuals as small bar charts overlaid on the network topology graph.
Each bus with measurements shows colored bars for each variable type (e.g., Vr, Vi, Pd, Qd),
with bar height proportional to the residual value.

# Arguments
- `SE_en`: State estimation results dictionary (with keys `"solution"`, `"objective"`, etc.)
- `math`: Mathematical model dictionary containing `"meas"`, `"bus"`, `"load"`, `"gen"`, etc.

# Keyword Arguments
- `phase::Int=1`: Terminal/phase to display (1, 2, or 3)
- `layout`: Graph layout function (default: `smart_layout`)
- `use_coords::Bool=false`: If `true`, use coordinate-based layout from bus lon/lat
- `bar_scale=:auto`: Scaling factor for bar heights (`:auto` computes from data)
- `bar_width=:auto`: Width of each bar (`:auto` computes from data)
- `show_node_labels::Bool=false`: Show bus ID labels
- `show_edge_labels::Bool=false`: Show edge labels
- `show_legend::Bool=true`: Show variable color legend
- `makie_backend`: Makie backend (default: `WGLMakie`)
- `figure_size::Tuple=(1000, 1200)`: Figure dimensions

# Returns
Tuple `(Figure, Axis, GraphPlot)`

# Example
```julia
fig, ax, gp = plot_residuals_bars(SE_result, math_model; phase=1, show_node_labels=true)
```
"""
function plot_residuals_bars(SE_en, math;
    phase::Int=1,
    layout=smart_layout,
    use_coords::Bool=false,
    bar_scale=:auto,
    bar_width=:auto,
    show_node_labels::Bool=false,
    show_edge_labels::Bool=false,
    show_legend::Bool=true,
    makie_backend=WGLMakie,
    figure_size=(1000, 1200),
    kwargs...)

    # Create network graph with appropriate layout
    if use_coords
        network_graph, coord_layout, _ = create_network_graph(math, layout)
        actual_layout = coord_layout
    else
        network_graph, _, _ = create_network_graph(math, layout)
        actual_layout = layout
    end

    # Collect residuals per bus
    bus_residuals = _collect_bus_residuals(SE_en, math)

    # Get node positions
    bus_positions = _get_node_positions(network_graph, actual_layout)

    # Plot base network
    fig, ax, gp = _base_network_plot(math, network_graph;
        layout=actual_layout, show_node_labels=show_node_labels,
        show_edge_labels=show_edge_labels,
        figure_size=figure_size, makie_backend=makie_backend, kwargs...)

    # Compute auto-scaling from position and residual ranges
    all_x = [p[1] for p in values(bus_positions)]
    all_y = [p[2] for p in values(bus_positions)]
    x_range = length(all_x) > 1 ? (maximum(all_x) - minimum(all_x)) : 1.0
    y_range = length(all_y) > 1 ? (maximum(all_y) - minimum(all_y)) : 1.0
    plot_scale = max(x_range, y_range)
    plot_scale = plot_scale == 0.0 ? 1.0 : plot_scale

    all_res_vals = Float64[]
    for (_, meas_list) in bus_residuals
        for m in meas_list
            val = get(m.res, phase, NaN)
            !isnan(val) && push!(all_res_vals, abs(val))
        end
    end
    max_res = isempty(all_res_vals) ? 1.0 : maximum(all_res_vals)
    max_res = max_res == 0.0 ? 1.0 : max_res

    actual_bar_scale = bar_scale === :auto ? (plot_scale * 0.15) / max_res : Float64(bar_scale)
    actual_bar_width = bar_width === :auto ? plot_scale * 0.015 : Float64(bar_width)

    seen_vars = Set{Symbol}()

    # Draw bars at each bus
    for (bus_id, meas_list) in bus_residuals
        !haskey(bus_positions, bus_id) && continue

        bx, by = bus_positions[bus_id]

        sorted_meas = sort(meas_list, by=m -> get(_VAR_ORDER, m.var, 99))
        n_vars = length(sorted_meas)

        total_width = n_vars * actual_bar_width * 1.5
        start_x = bx - total_width / 2 + actual_bar_width * 0.75

        for (i, m) in enumerate(sorted_meas)
            val = get(m.res, phase, NaN)
            isnan(val) && continue

            push!(seen_vars, m.var)

            bar_x = start_x + (i - 1) * actual_bar_width * 1.5
            bar_h = abs(val) * actual_bar_scale
            bar_color = get(_VAR_COLORS, m.var, :gray)

            offset_y = plot_scale * 0.02
            rect_points = Point2f[
                (bar_x - actual_bar_width/2, by + offset_y),
                (bar_x + actual_bar_width/2, by + offset_y),
                (bar_x + actual_bar_width/2, by + offset_y + bar_h),
                (bar_x - actual_bar_width/2, by + offset_y + bar_h),
            ]
            poly!(ax, rect_points; color=bar_color, strokecolor=:black, strokewidth=0.5)
        end
    end

    # Add legend
    if show_legend && !isempty(seen_vars)
        sorted_vars = sort(collect(seen_vars), by=v -> get(_VAR_ORDER, v, 99))
        legend_entries = [PolyElement(color=get(_VAR_COLORS, v, :gray), strokecolor=:black, strokewidth=0.5) for v in sorted_vars]
        legend_labels = [get(_VAR_LABELS, v, string(v)) for v in sorted_vars]
        Legend(fig[1, 2], legend_entries, legend_labels, "Variables")
    end

    return fig, ax, gp
end


"""
    plot_residuals_bubbles(SE_en, math; kwargs...)

Plot measurement residuals as circles overlaid on the network topology graph.
Circle radius is proportional to the aggregated residual magnitude at each bus,
and color is mapped via a continuous colormap.

# Arguments
- `SE_en`: State estimation results dictionary
- `math`: Mathematical model dictionary with measurements

# Keyword Arguments
- `phase::Int=1`: Terminal/phase to display (1, 2, or 3)
- `layout`: Graph layout function (default: `smart_layout`)
- `use_coords::Bool=false`: If `true`, use coordinate-based layout
- `max_bubble_size::Float64=30.0`: Maximum bubble marker size in pixels
- `min_bubble_size::Float64=5.0`: Minimum bubble marker size in pixels
- `aggregation::Symbol=:max`: How to aggregate residuals per bus (`:max`, `:mean`, `:sum`)
- `show_node_labels::Bool=false`: Show bus ID labels
- `show_legend::Bool=true`: Show colorbar
- `colormap=:RdYlGn_r`: Colormap for bubble coloring
- `makie_backend`: Makie backend (default: `WGLMakie`)
- `figure_size::Tuple=(1000, 1200)`: Figure dimensions

# Returns
Tuple `(Figure, Axis, GraphPlot)`

# Example
```julia
fig, ax, gp = plot_residuals_bubbles(SE_result, math_model; aggregation=:mean, colormap=:viridis)
```
"""
function plot_residuals_bubbles(SE_en, math;
    phase::Int=1,
    layout=smart_layout,
    use_coords::Bool=false,
    max_bubble_size::Float64=30.0,
    min_bubble_size::Float64=5.0,
    aggregation::Symbol=:max,
    show_node_labels::Bool=false,
    show_edge_labels::Bool=false,
    show_legend::Bool=true,
    colormap=:RdYlGn_r,
    makie_backend=WGLMakie,
    figure_size=(1000, 1200),
    kwargs...)

    # Create network graph
    if use_coords
        network_graph, coord_layout, _ = create_network_graph(math, layout)
        actual_layout = coord_layout
    else
        network_graph, _, _ = create_network_graph(math, layout)
        actual_layout = layout
    end

    # Collect and aggregate residuals
    bus_residuals = _collect_bus_residuals(SE_en, math)
    bus_positions = _get_node_positions(network_graph, actual_layout)
    bus_agg = _aggregate_bus_residuals(bus_residuals, phase, aggregation)

    # Plot base network
    fig, ax, gp = _base_network_plot(math, network_graph;
        layout=actual_layout, show_node_labels=show_node_labels,
        show_edge_labels=show_edge_labels,
        figure_size=figure_size, makie_backend=makie_backend, kwargs...)

    if !isempty(bus_agg)
        agg_values = collect(values(bus_agg))
        min_val = minimum(agg_values)
        max_val = maximum(agg_values)
        val_range = max_val - min_val
        val_range = val_range == 0.0 ? 1.0 : val_range

        bubble_x = Float64[]
        bubble_y = Float64[]
        bubble_sizes = Float64[]
        bubble_colors = Float64[]

        for (bus_id, agg_val) in bus_agg
            !haskey(bus_positions, bus_id) && continue
            bx, by = bus_positions[bus_id]
            push!(bubble_x, bx)
            push!(bubble_y, by)

            norm_val = (agg_val - min_val) / val_range
            push!(bubble_sizes, min_bubble_size + norm_val * (max_bubble_size - min_bubble_size))
            push!(bubble_colors, agg_val)
        end

        scatter!(ax, bubble_x, bubble_y;
            markersize=bubble_sizes,
            color=bubble_colors,
            colormap=colormap,
            colorrange=(min_val, max_val),
            strokecolor=:black,
            strokewidth=0.5,
            alpha=0.7)

        if show_legend
            Colorbar(fig[1, 2]; colormap=colormap, colorrange=(min_val, max_val),
                label="Residual ($(aggregation))")
        end
    end

    return fig, ax, gp
end


"""
    plot_residuals_heatmap(SE_en, math; kwargs...)

Plot the network graph with nodes colored by their aggregated measurement
residual values using a continuous colormap. Nodes without measurements
are shown in translucent gray.

# Arguments
- `SE_en`: State estimation results dictionary
- `math`: Mathematical model dictionary with measurements

# Keyword Arguments
- `phase::Int=1`: Terminal/phase to display (1, 2, or 3)
- `layout`: Graph layout function (default: `smart_layout`)
- `use_coords::Bool=false`: If `true`, use coordinate-based layout
- `aggregation::Symbol=:max`: How to aggregate residuals per bus (`:max`, `:mean`, `:sum`)
- `show_node_labels::Bool=false`: Show bus ID labels
- `show_legend::Bool=true`: Show colorbar
- `colormap=:RdYlGn_r`: Colormap for node coloring
- `node_size::Int=15`: Size of nodes that have measurements
- `makie_backend`: Makie backend (default: `WGLMakie`)
- `figure_size::Tuple=(1000, 1200)`: Figure dimensions

# Returns
Tuple `(Figure, Axis, GraphPlot)`

# Example
```julia
fig, ax, gp = plot_residuals_heatmap(SE_result, math_model; colormap=:hot, aggregation=:mean)
```
"""
function plot_residuals_heatmap(SE_en, math;
    phase::Int=1,
    layout=smart_layout,
    use_coords::Bool=false,
    aggregation::Symbol=:max,
    show_node_labels::Bool=false,
    show_edge_labels::Bool=false,
    show_legend::Bool=true,
    colormap=:RdYlGn_r,
    node_size::Int=15,
    makie_backend=WGLMakie,
    figure_size=(1000, 1200),
    kwargs...)

    # Create network graph
    if use_coords
        network_graph, coord_layout, _ = create_network_graph(math, layout)
        actual_layout = coord_layout
    else
        network_graph, _, _ = create_network_graph(math, layout)
        actual_layout = layout
    end

    # Collect and aggregate residuals
    bus_residuals = _collect_bus_residuals(SE_en, math)
    bus_agg = _aggregate_bus_residuals(bus_residuals, phase, aggregation)

    # Compute node colors based on residuals
    agg_values = isempty(bus_agg) ? [0.0] : collect(values(bus_agg))
    min_val = minimum(agg_values)
    max_val = maximum(agg_values)
    val_range = max_val - min_val
    val_range = val_range == 0.0 ? 1.0 : val_range

    cmap = Makie.to_colormap(colormap)
    n_colors = length(cmap)

    node_colors = []
    node_sizes = Int[]
    for v in 1:nv(network_graph)
        bus_id = string(get_prop(network_graph, v, :bus_id))
        if haskey(bus_agg, bus_id)
            norm_val = (bus_agg[bus_id] - min_val) / val_range
            color_idx = clamp(round(Int, norm_val * (n_colors - 1)) + 1, 1, n_colors)
            push!(node_colors, cmap[color_idx])
            push!(node_sizes, node_size)
        else
            push!(node_colors, Makie.RGBAf(0.7, 0.7, 0.7, 0.5))
            push!(node_sizes, 5)
        end
    end

    # Plot with overridden node colors
    fig, ax, gp = _base_network_plot(math, network_graph;
        layout=actual_layout, show_node_labels=show_node_labels,
        show_edge_labels=show_edge_labels,
        figure_size=figure_size, makie_backend=makie_backend,
        node_color_override=node_colors,
        node_size_override=node_sizes,
        kwargs...)

    if show_legend && !isempty(bus_agg)
        Colorbar(fig[1, 2]; colormap=colormap, colorrange=(min_val, max_val),
            label="Residual ($(aggregation))")
    end

    return fig, ax, gp
end


# Re-export PMDSE utility functions from parent module
export viz_residuals
export df_meas_res
export add_pd_qd_vmn!
export write_sm_measurements
export plot_residuals_bars
export plot_residuals_bubbles
export plot_residuals_heatmap









end # module PMDSEUtils


