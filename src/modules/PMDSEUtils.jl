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

    MAPE, APEs_df, _ = Pliers._calculate_MAPE(SE_en, PF_en, math_en)
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
    extra_keys(meass, ["cmp", "cmp_id", "var", "crit", "res"])
        
end


function write_sm_measurements(PF_RES, math, measurements_file; σ=0.05, measurement_model =PowerModelsDistributionStateEstimation.IndustrialENMeasurementsModel )
    dictify_solution!(PF_RES["solution"], math)
    math_meas_en = add_vmn_p_q(math, PF_RES["solution"])
    write_measurements!(measurement_model, math_meas_en, PF_RES, measurements_file, σ=σ) 
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

















# Re-export PMDSE utility functions from parent module
export viz_residuals
export df_meas_res
export add_pd_qd_vmn!
export write_sm_measurements









end # module PMDSEUtils


