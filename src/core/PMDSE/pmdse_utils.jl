
function viz_residuals(SE_en, math_en, se_sol_en)
        
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

    min_val = minimum([
        minimum(filter(!isnan, df_meas_res.Res1)), 
        minimum(filter(!isnan, df_meas_res.Res2)), 
        minimum(filter(!isnan, df_meas_res.Res3))
    ])

    max_val = maximum([
        maximum(filter(!isnan, df_meas_res.Res1)), 
        maximum(filter(!isnan, df_meas_res.Res2)), 
        maximum(filter(!isnan, df_meas_res.Res3))
    ])

    # Define color intervals
    num_intervals = 6

    lowerBounds = [[min_val + i * (0.01 - min_val) / 3 for i in 0:3]..., 1e-1, 1]
    uppoerBounds = [[0.00333332833981672 + i * (0.1 - 0.00333332833981672) / 3 for i in 0:3]..., 1, max_val+1e-6]


    colors = [crayon"green bold", crayon"cyan bold", crayon"blue bold", crayon"magenta bold", crayon"yellow bold", crayon"red bold"]


    # Create highlighters for each interval
    highlighters = []
    for i in 1:num_intervals
        lower_bound = lowerBounds[i]
        upper_bound = uppoerBounds[i]
        push!(highlighters, Highlighter(
            (data, row, col) -> col in 2:4 && data[row, col] >= lower_bound && data[row, col] < upper_bound,
            colors[i]
        ))
    end
    hideNan_hl = Highlighter(       
        (data, i, j) -> (j ∈ collect(2:4) && isnan(data[i,j])),
        crayon"dark_gray conceal"
        );
    push!(highlighters, hideNan_hl);

    # Define headers
    header = ["Meas", "Res1", "Res2", "Res3"]
    res_heatmap = pretty_table(df_meas_res, header=header, header_crayon=crayon"fg:yellow", highlighters=Tuple(highlighters), tf=tf_unicode)

    df_legend = DataFrame(
        Lower_Bound = lowerBounds,
        Upper_Bound = uppoerBounds,
    )

    legend_highlighters = []
    for i in 1:num_intervals
        push!(legend_highlighters, Highlighter(
            (data, row, col) -> (row == i),
            colors[i]
        ))
    end
    println("Legend:")
    legend = pretty_table(df_legend, header=["Lower Bound", "Upper Bound"], header_crayon=crayon"fg:yellow", highlighters= Tuple(legend_highlighters), tf=tf_unicode)

    display(res_heatmap)
    println(legend)

end

#TODO: fix it to not only show 3 columns but extend the res array,, also note math has maybe all things u need?!
function viz_residuals(SE_en, math_en)

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

    min_val = minimum([
        minimum(filter(!isnan, df_meas_res.Res1)), 
        minimum(filter(!isnan, df_meas_res.Res2)), 
        minimum(filter(!isnan, df_meas_res.Res3))
    ])

    max_val = maximum([
        maximum(filter(!isnan, df_meas_res.Res1)), 
        maximum(filter(!isnan, df_meas_res.Res2)), 
        maximum(filter(!isnan, df_meas_res.Res3))
    ])

    # Define color intervals
    num_intervals = 6

    lowerBounds = [[min_val + i * (0.01 - min_val) / 3 for i in 0:3]..., 1e-1, 1]
    uppoerBounds = [[0.00333332833981672 + i * (0.1 - 0.00333332833981672) / 3 for i in 0:3]..., 1, max_val+1e-6]


    colors = [crayon"green bold", crayon"cyan bold", crayon"blue bold", crayon"magenta bold", crayon"yellow bold", crayon"red bold"]


    # Create highlighters for each interval
    highlighters = []
    for i in 1:num_intervals
        lower_bound = lowerBounds[i]
        upper_bound = uppoerBounds[i]
        push!(highlighters, Highlighter(
            (data, row, col) -> col in 2:4 && data[row, col] >= lower_bound && data[row, col] < upper_bound,
            colors[i]
        ))
    end
    hideNan_hl = Highlighter(       
        (data, i, j) -> (j ∈ collect(2:4) && isnan(data[i,j])),
        crayon"dark_gray conceal"
        );
    push!(highlighters, hideNan_hl);

    # Define headers
    header = ["Meas", "Res1", "Res2", "Res3"]
    res_heatmap = pretty_table(df_meas_res, header=header, header_crayon=crayon"fg:yellow", highlighters=Tuple(highlighters), tf=tf_unicode)

    df_legend = DataFrame(
        Lower_Bound = lowerBounds,
        Upper_Bound = uppoerBounds,
    )

    legend_highlighters = []
    for i in 1:num_intervals
        push!(legend_highlighters, Highlighter(
            (data, row, col) -> (row == i),
            colors[i]
        ))
    end
    println("Legend:")
    legend = pretty_table(df_legend, header=["Lower Bound", "Upper Bound"], header_crayon=crayon"fg:yellow", highlighters= Tuple(legend_highlighters), tf=tf_unicode)

    display(res_heatmap)
    println(legend)

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

    se_sol  = SE_RES["solution"]
    pf_sol = PF_RES["solution"]

    dictify_solution!(pf_sol, math)
    dictify_solution!(se_sol, math)

    APEs = [] 
    Errors = [] 

    Errors_df = DataFrame(Bus = String[], Terminal = String[], Error = ComplexF64[])
    APEs_df = DataFrame(Bus = String[], Terminal = String[], APE = Float64[])
    for (b,bus) in se_sol["bus"]

        for (term,Vse) in bus["voltage"]
            # println("Bus    :",b)
            # println("Term   :",term)
            # println("Vse   :",Vse)
            # println("Vpf     :", pf_sol["bus"][b]["voltage"][term])
            Vpf = pf_sol["bus"][b]["voltage"][term]

            # APE = abs.(Vpf) == 0 ? 0 : abs.( Vse - Vpf ) ./ abs.(Vpf) * 100
            APE = ( abs.( Vse - Vpf ) ./ abs.(Vpf) )* 100
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
