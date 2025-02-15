
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

