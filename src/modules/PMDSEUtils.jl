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
using ..Pliers.PMDGraph: Graphs, MetaGraphs, create_network_graph, smart_layout,
    network_graph_plot, _decorate_nodes!, _decorate_edges!
using Graphs: vertices, nv, edges
using MetaGraphs: get_prop, props, set_prop!

# plotting packages
using ..Pliers: Makie, CairoMakie, WGLMakie
using Makie: Point2f, poly!, PolyElement, Legend, Colorbar, scatter!

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
function viz_residuals(SE_en, math_en; detailed = true, show_legend = false, mape = nothing, mape_vmn = nothing, APEs_df = nothing)#pm_form = PowerModelsDistribution.IVRENPowerModel)

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
    !isnothing(mape)     ? (isapprox(mape,     0, atol=0.05) ? printstyled(" MAPE_Vg  (phase-to-gnd) : $(round(mape,     digits=4)) %\n", color=:green) : printstyled(" MAPE_Vg  (phase-to-gnd) : $(round(mape,     digits=4)) %\n", color=:red)) : nothing
    !isnothing(mape_vmn) ? (isapprox(mape_vmn, 0, atol=0.05) ? printstyled(" MAPE_Vmn (phase-to-ntr) : $(round(mape_vmn, digits=4)) %\n", color=:green) : printstyled(" MAPE_Vmn (phase-to-ntr) : $(round(mape_vmn, digits=4)) %\n", color=:red)) : nothing
    
    m,n, dof = get_m_n_dof(math_en)
    printstyled(" m : $(m) \n", color=:green)
    printstyled(" n : $(n) \n", color=:green)
    dof > 0 ? printstyled(" Degrees of freedom : $(dof) \n", color=:green) : printstyled(" Degrees of freedom : $(dof) \n", color=:red)
    
    # @suppress display(instantiate_mc_model(
    #     math_en,
    #     pm_form,
    #     PowerModelsDistributionStateEstimation.build_mc_se;
    # ).model)

        if detailed
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
        header_names = ["Meas", "Res1", "Res2", "Res3", "Res4"]
        try
            res_heatmap = pretty_table(df_meas_res; column_labels=header_names, highlighters=highlighters)
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
    end

    return solve_time, objective, termination_status, primal_status, mape, mape_vmn, m, n, dof

end

function viz_residuals(SE_en, math_en, PF_en; detailed=true, show_legend=false, kwargs...)

    mape     = _calculate_MAPE(SE_en, PF_en, math_en)
    mape_vmn = _calculate_MAPE_vmn(SE_en, PF_en, math_en)
    viz_residuals(SE_en, math_en; detailed = detailed, show_legend = show_legend, mape = mape, mape_vmn = mape_vmn, kwargs...)

end

"""
    _get_bus_voltages(sol_dict, math) → Dict{String, Dict{String, ComplexF64}}

Extract bus voltage phasors keyed by bus_id → terminal_string → complex voltage.
Handles both raw `solution` dicts and pre-dictified solution dicts.
"""
function _get_bus_voltages(sol_dict::Dict{String,Any}, math::Dict{String,Any})
    buses = haskey(sol_dict, "bus") ? sol_dict["bus"] : sol_dict["solution"]["bus"]
    result = Dict{String, Dict{String, ComplexF64}}()
    for (b, bus_data) in buses
        if haskey(bus_data, "voltage")
            result[b] = Dict{String, ComplexF64}(
                string(k) => ComplexF64(v) for (k, v) in bus_data["voltage"]
            )
        elseif haskey(bus_data, "vr") && haskey(bus_data, "vi")
            # raw vr/vi arrays — terminals inferred from math
            terms = math["bus"][b]["terminals"]
            vr = bus_data["vr"]
            vi = bus_data["vi"]
            volt_dict = Dict{String, ComplexF64}()
            for (i, t) in enumerate(terms)
                if i <= length(vr)
                    volt_dict[string(t)] = ComplexF64(vr[i], vi[i])
                end
            end
            result[b] = volt_dict
        end
    end
    return result
end

"""
    _calculate_MAPE(se_sol, pf_sol, math) → Float64

Mean Absolute Percentage Error for **phase-to-ground** (phase-to-zero-reference)
voltage magnitudes, phases 1–3 only. Terminal 4 (neutral) is excluded from
the comparison set so that EN and KRN results are directly comparable.

    MAPE_Vg = mean over all (bus, phase) of ||V_SE| - |V_PF|| / |V_PF| × 100 %

For the phase-to-neutral analogue see `_calculate_MAPE_vmn`.

Both `se_sol` and `pf_sol` accept either a raw result dict (with a `"solution"`
key) or a pre-dictified solution dict.
"""
function _calculate_MAPE(
    se_sol::Dict{String,Any},
    pf_sol::Dict{String,Any},
    math::Dict{String,Any},
)::Float64

    se_voltages = _get_bus_voltages(se_sol, math)
    pf_voltages = _get_bus_voltages(pf_sol, math)

    apes = Float64[]
    for (b, se_volts) in se_voltages
        haskey(pf_voltages, b) || continue
        pf_volts = pf_voltages[b]
        for phase in ("1", "2", "3")
            (haskey(se_volts, phase) && haskey(pf_volts, phase)) || continue
            Vse = abs(se_volts[phase])
            Vpf = abs(pf_volts[phase])
            Vpf == 0.0 && continue
            push!(apes, abs(Vse - Vpf) / Vpf * 100.0)
        end
    end
    isempty(apes) && return NaN
    return mean(filter(!isnan, apes))
end

"""
    _calculate_MAPE_vmn(se_sol, pf_sol, math) → Float64

Mean Absolute Percentage Error for **phase-to-neutral** voltage magnitudes,
phases 1–3 only (terminal 4 excluded from the comparison set).

    MAPE_Vmn = mean over all (bus, phase) of ||V_SE_mn| - |V_PF_mn|| / |V_PF_mn| × 100 %

where ``V_{mn} = V_{\\text{phase}} - V_{\\text{neutral}}``.

Neutral reference handling:
  - **EN models** (terminal `"4"` present in bus voltage dict): `V_neutral = V["4"]`
  - **KRN models** (terminal `"4"` absent): `V_neutral = 0 + 0im`
    (neutral is the network reference node; V_phase already equals V_phase-to-neutral)

This convention means KRN SE results are correctly compared against EN PF
ground truth: the neutral displacement `|V_4_pf|` is subtracted from the PF
phasors before the APE is computed, so only the estimation error contributes.

Both `se_sol` and `pf_sol` accept either a raw result dict (with a `"solution"`
key) or a pre-dictified solution dict.
"""
function _calculate_MAPE_vmn(
    se_sol::Dict{String,Any},
    pf_sol::Dict{String,Any},
    math::Dict{String,Any},
)::Float64

    se_voltages = _get_bus_voltages(se_sol, math)
    pf_voltages = _get_bus_voltages(pf_sol, math)

    apes = Float64[]
    for (b, se_volts) in se_voltages
        haskey(pf_voltages, b) || continue
        pf_volts = pf_voltages[b]
        # Neutral voltage: 0+0im when absent (Kron: neutral is reference = 0)
        V4_se = get(se_volts, "4", ComplexF64(0))
        V4_pf = get(pf_volts, "4", ComplexF64(0))
        for phase in ("1", "2", "3")
            (haskey(se_volts, phase) && haskey(pf_volts, phase)) || continue
            Vse_mn = abs(se_volts[phase] - V4_se)
            Vpf_mn = abs(pf_volts[phase] - V4_pf)
            Vpf_mn == 0.0 && continue
            push!(apes, abs(Vse_mn - Vpf_mn) / Vpf_mn * 100.0)
        end
    end
    isempty(apes) && return NaN
    return mean(filter(!isnan, apes))
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
        # Handle cases where layout_fn returns Point2f/Point2 instead of Tuple
        pos = positions[v]
        bus_positions[bus_id] = (Float64(pos[1]), Float64(pos[2]))
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
                Statistics.mean(vals)
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
- `colormap=Makie.Reverse(:RdYlGn)`: Colormap for bubble coloring
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
    colormap=Makie.Reverse(:RdYlGn),
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
- `colormap=Makie.Reverse(:RdYlGn)`: Colormap for node coloring
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
    colormap=Makie.Reverse(:RdYlGn),
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


# ═══════════════════════════════════════════════════════════════════════════════
# Measurement Overlay Visualization on Network Graph
# ═══════════════════════════════════════════════════════════════════════════════

# Classify each measurement variable into a semantic group
const _VAR_GROUP = Dict{Symbol, Symbol}(
    :vmn => :voltage, :vm  => :voltage, :vr  => :voltage, :vi  => :voltage,
    :vll => :voltage, :w   => :voltage,
    :pd  => :power,   :qd  => :power,   :pg  => :power,   :qg  => :power,
    :ptot => :power,  :qtot => :power, :pd_bus => :power, :qd_bus => :power,
    :crd => :current, :cid => :current, :crg => :current, :cig => :current,
    :crd_bus => :current, :cid_bus => :current,
)

"""
Node color based on which measurement groups (V / P-Q / I) are present at bus.

| Color       | Groups present           |
|-------------|--------------------------|
| deepskyblue | voltage only             |
| forestgreen | power only               |
| crimson     | current only             |
| darkorchid  | voltage + power          |
| teal        | voltage + current        |
| saddlebrown | power + current          |
| gold        | voltage + power + current|
"""
function _coverage_color(groups::Set{Symbol})
    has_v = :voltage ∈ groups
    has_p = :power   ∈ groups
    has_i = :current ∈ groups
    has_v && has_p && has_i && return :gold
    has_v && has_p            && return :darkorchid
    has_v && has_i            && return :teal
    has_p && has_i            && return :saddlebrown
    has_v                     && return :deepskyblue
    has_p                     && return :forestgreen
    has_i                     && return :crimson
    return :gray
end

function _coverage_label(groups::Set{Symbol})
    parts = String[]
    :voltage ∈ groups && push!(parts, "V")
    :power   ∈ groups && push!(parts, "P/Q")
    :current ∈ groups && push!(parts, "I")
    return isempty(parts) ? "other" : join(parts, "+")
end

"""Resolve which bus a measurement is attached to (returns bus_id string or nothing)."""
function _meas_bus_id(meas::Dict, math::Dict)
    cmp    = meas["cmp"]
    cmp_id = string(meas["cmp_id"])
    cmp == :bus  && return cmp_id
    cmp == :load && return string(math["load"][cmp_id]["load_bus"])
    cmp == :gen  && return string(math["gen"][cmp_id]["gen_bus"])
    return nothing
end

"""Compact single-line description of one measurement (fits tooltip line ≤60 chars)."""
function _meas_compact_str(meas_id::String, meas::Dict)
    var_str = string(meas["var"])
    cmp_str = "$(meas["cmp"])[$(meas["cmp_id"])]"
    dsts    = meas["dst"]
    dst_str = try
        μs = [round(d.μ, digits=4) for d in dsts]
        σ1 = round(dsts[1].σ, digits=4)
        length(dsts) == 1 ? "μ=$(μs[1]) σ=$σ1" : "μ=$(μs) σ=$σ1"
    catch
        join(string.(dsts), " | ")
    end
    res_str = (haskey(meas, "res") && !isnothing(meas["res"])) ?
        " res=$(round.(meas["res"], digits=3))" : ""
    return "[$meas_id] $var_str@$cmp_str  $dst_str$res_str"
end


"""
    plot_measurements(math; kwargs...)

Overlay measurement indicators from `math["meas"]` directly on the network
topology graph.  Each bus that carries at least one measurement is re-drawn
with a **hexagon** marker (meter symbol) whose color encodes which measurement
groups are present:

| Color        | Coverage              |
|--------------|-----------------------|
| deepskyblue  | V only                |
| forestgreen  | P/Q only              |
| crimson      | I only                |
| darkorchid   | V + P/Q               |
| teal         | V + I                 |
| saddlebrown  | P/Q + I               |
| gold         | V + P/Q + I           |

Hovering over a measured node shows a DataInspector tooltip with one compact
line per measurement, and simultaneously prints the full measurement details
to the terminal.

# Arguments
- `math`: MATH-format dict with a `"meas"` key.

# Keyword Arguments
- `use_coords::Bool=false`: geographic coordinate layout
- `layout`: graph layout function (default: `smart_layout`)
- `meter_size::Int=22`: marker size for measured nodes (px)
- `show_node_labels::Bool=false`: show bus ID labels
- `show_edge_labels::Bool=false`: show edge labels
- `show_legend::Bool=true`: add coverage legend
- `makie_backend`: Makie backend (default: `WGLMakie`)
- `figure_size::Tuple=(1000,1200)`: figure dimensions

# Returns
`(Figure, Axis, GraphPlot)`
"""
function plot_measurements(math;
    use_coords::Bool=false,
    layout=smart_layout,
    meter_size::Int=22,
    show_node_labels::Bool=false,
    show_edge_labels::Bool=false,
    show_legend::Bool=true,
    makie_backend=WGLMakie,
    figure_size=(1000, 1200),
    kwargs...)

    haskey(math, "meas") || error("No `meas` key in math dict")

    # Build graph
    if use_coords
        network_graph, coord_layout, _ = create_network_graph(math, layout)
        actual_layout = coord_layout
    else
        network_graph, _, _ = create_network_graph(math, layout)
        actual_layout = layout
    end

    # Collect measurements per bus (sorted by meas id)
    bus_meas = Dict{String, Vector{Tuple{String,Dict}}}()
    for (meas_id, meas) in sort(collect(math["meas"]),
            by = x -> something(tryparse(Int, x[1]), typemax(Int)))
        bus_id = _meas_bus_id(meas, math)
        isnothing(bus_id) && continue
        push!(get!(bus_meas, bus_id, Tuple{String,Dict}[]), (meas_id, meas))
    end

    # Apply default decoration, then override for measured buses
    _decorate_nodes!(network_graph, math)
    _decorate_edges!(network_graph, math)

    seen_coverages = Dict{String, Symbol}()  # coverage_label → color (for legend)

    for v in vertices(network_graph)
        bus_id = string(get_prop(network_graph, v, :bus_id))
        !haskey(bus_meas, bus_id) && continue

        meass  = bus_meas[bus_id]
        groups = Set{Symbol}(get(_VAR_GROUP, meas["var"], :unknown) for (_, meas) in meass)
        color  = _coverage_color(groups)
        label  = _coverage_label(groups)
        seen_coverages[label] = color

        # Override visual properties
        vp = network_graph.vprops[v]
        vp[:node_color]  = color
        vp[:node_marker] = :hexagon
        vp[:marker_size] = meter_size

        # Store measurement summaries as vertex properties for DataInspector
        vp[:meas_coverage] = label
        vp[:meas_count]    = length(meass)
        for (meas_id, meas) in meass
            vp[Symbol("meas_$meas_id")] = _meas_compact_str(meas_id, meas)
        end
    end

    # Extract per-node/edge arrays from the now-decorated graph
    node_color   = [get_prop(network_graph, i, :node_color)  for i in 1:nv(network_graph)]
    node_marker  = [get_prop(network_graph, i, :node_marker) for i in 1:nv(network_graph)]
    node_size    = [get_prop(network_graph, i, :marker_size) for i in 1:nv(network_graph)]
    edge_color   = [get_prop(network_graph, e, :edge_color)   for e in edges(network_graph)]
    arrow_marker = [get_prop(network_graph, e, :arrow_marker) for e in edges(network_graph)]
    arrow_size   = [get_prop(network_graph, e, :arrow_size)   for e in edges(network_graph)]
    arrow_shift  = [get_prop(network_graph, e, :arrow_shift)  for e in edges(network_graph)]

    nlabels = show_node_labels ?
        [string(get_prop(network_graph, i, :bus_id)) for i in 1:nv(network_graph)] : nothing

    # Plot via network_graph_plot — DataInspector is activated inside via
    # _build_inspector_labels, which reads all vertex props (incl. meas_*) for tooltip
    fig, ax, gp = network_graph_plot(
        network_graph;
        layout=actual_layout,
        figure_size=figure_size,
        makie_backend=makie_backend,
        show_node_labels=show_node_labels,
        nlabels=nlabels,
        show_edge_labels=show_edge_labels,
        node_color=node_color,
        node_marker=node_marker,
        node_size=node_size,
        edge_color=edge_color,
        arrow_show=true,
        arrow_marker=arrow_marker,
        arrow_size=arrow_size,
        arrow_shift=arrow_shift,
        _pmd_data=math,
        kwargs...)

    # Coverage legend
    if show_legend && !isempty(seen_coverages)
        sorted = sort(collect(seen_coverages), by=x -> x[1])
        Legend(fig[1, 2],
            [PolyElement(color=c, strokecolor=:black, strokewidth=0.5) for (_, c) in sorted],
            [lbl for (lbl, _) in sorted],
            "Measurements")
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
export plot_measurements









end # module PMDSEUtils
