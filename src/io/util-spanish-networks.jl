"""
    rm_spanish_transformer!(data_eng)

Remove the Spanish transformer from the given `data_eng` dictionary.

# Arguments
- `data_eng`: A dictionary containing the network data.

# Description
This function removes the Spanish transformer from the `data_eng` dictionary. It updates the voltage source parameters, deletes the transformer, deletes the lines connected before the transformer, updates the bus settings, and resets the terminals and grounded status of the buses.

"""
function rm_spanish_transformer!(data_eng)
   
    if haskey(data_eng, "transformer")
        trans = first(data_eng["transformer"]).second

        old_slack = first(data_eng["transformer"]).second["bus"][1]
        new_slack = first(data_eng["transformer"]).second["bus"][2]
        #display("Old slack bus: $old_slack")
        #display("New slack bus: $new_slack")

        vprim_scale = trans["vm_nom"][2]/trans["vm_nom"][1]
        data_eng["voltage_source"]["source"]["vm"] *= vprim_scale
        data_eng["voltage_source"]["source"]["rs"] *= vprim_scale^2
        data_eng["voltage_source"]["source"]["xs"] *= vprim_scale^2
        data_eng["voltage_source"]["source"]["bus"] =   new_slack # the new slack bus

        delete!(data_eng, "transformer")
        
        #delete lines connected before transformer
        
        line_bf = string.([x for(x, line) in data_eng["line"] if line["f_bus"] == old_slack || line["t_bus"] == old_slack]...)
        delete!(data_eng["line"], line_bf)

        delete!(data_eng["bus"], old_slack)
        data_eng["settings"]["vbases_default"][new_slack] = data_eng["settings"]["vbases_default"]["source"]*vprim_scale
        #delete!(data_eng["settings"]["vbases_default"], "source")
        delete!(data_eng["settings"]["vbases_default"], "source")
        delete!(data_eng["bus"], "source")

        for (b, bus) in data_eng["bus"]
            if length(bus["terminals"]) > 4
                bus["terminals"] = [1, 2, 3, 4]
                bus["grounded"] = Int64[]
            end
        end

        for (sh, shunt) in data_eng["shunt"]
            #if shunt["bus"] == new_slack
                delete!(data_eng["shunt"], sh)
            #end
        end

    end
end




"""
    add_needed_details!(data_eng)

Add additional details to the `data_eng` dictionary for lines and loads.

# Arguments
- `data_eng`: A dictionary containing engineering data.

# Description
This function adds the following details to the `data_eng` dictionary for lines:
- `"Rs (Ω)"`: The resistance of the line calculated as the product of the line length and the resistance per unit length (`rs`).
- `"Xs (Ω)"`: The reactance of the line calculated as the product of the line length and the reactance per unit length (`xs`).

For loads, the function adds the following details:
- `"PF"`: The power factor calculated as the ratio of the active power (`pd_nom`) to the apparent power (`sqrt(pd_nom^2 + qd_nom^2)`).
- `"S (VAh)"`: The apparent power calculated as the square root of the sum of the squares of the active power and reactive power (`sqrt(pd_nom^2 + qd_nom^2)`).
- `"connected_phase"`: The phase connection of the load, represented as a string (`"a"`, `"b"`, `"c"`, or `"G"`) based on the value of the `connections` field.

"""
function add_needed_details!(data_eng)
    if haskey(data_eng, "line")
        
        for (l, line) in data_eng["line"]
            linecode = data_eng["linecode"][line["linecode"]]
            line["Rs (Ω)"] =  line["length"]*linecode["rs"]
            line["Xs (Ω)"] =  line["length"]*linecode["xs"]
        end

        for (lo, load) in data_eng["load"]
            load["PF"] = load["pd_nom"] / sqrt((load["pd_nom"][1])^2 + (load["qd_nom"][1])^2)
            load["S (VAh)"] =  sqrt((load["pd_nom"][1])^2 + (load["qd_nom"][1])^2)
            load["connected_phase"] = ["a" "b" "c" "G"][load["connections"][1]]
        end
    end
end


function add_coordinates!(data_eng, locs_x, locs_y)
    if haskey(data_eng, "bus")
        for (b, bus) in data_eng["bus"]
            
            if b != "sourcebus"
                b = parse(Int64, b)
                bus["X_coord"] = locs_x[b]
                bus["Y_coord"] = locs_y[b]
            else
                bus["X_coord"] = maximum(locs_x)
                bus["Y_coord"] = minimum(locs_y)
            end
        end
    end
    
end

"""
    attach_loads!(data_eng)

Attach loads to the specified bus in the given `data_eng` dictionary.

# Arguments
- `data_eng`: A dictionary containing network data.

"""
function attach_loads!(data_eng)
    if haskey(data_eng, "load")
        for (l, loadelement) in data_eng["load"]
            if haskey(loadelement, "bus")
                bus_idx = loadelement["bus"]
                data_eng["bus"][bus_idx]["load_S"] = loadelement["S (VAh)"]
                data_eng["bus"][bus_idx]["load_PF"] = loadelement["PF"]
                data_eng["bus"][bus_idx]["load_connected_phase"] = loadelement["connected_phase"]
                data_eng["bus"][bus_idx]["load_pd_nom"] = loadelement["pd_nom"]
                data_eng["bus"][bus_idx]["load_qd_nom"] = loadelement["qd_nom"]
                data_eng["bus"][bus_idx]["configuration"] = loadelement["configuration"]
                data_eng["bus"][bus_idx]["loaded"] = true
            end
        end
    end
end

"""
    attach_linecodes!(data_eng)

Attach linecodes to the given network data.

# Arguments
- `data_eng`: A dictionary containing network data.

# Description
This function attaches linecodes to the lines in the network data. It retrieves the linecode information from the `linecode` dictionary and assigns the corresponding values to the line dictionary.

"""
function attach_linecodes!(data_eng)
    if haskey(data_eng, "line")
        for (l, line) in data_eng["line"]
            linecode = data_eng["linecode"][line["linecode"]]
            line["lc_b_fr"] = linecode["b_fr"]
            line["lc_rs"] = linecode["rs"]
            line["lc_cm_ub"] = linecode["cm_ub"]
            line["lc_xs"] = linecode["xs"]
            line["lc_b_to"] = linecode["b_to"]
            line["lc_g_to"] = linecode["g_to"]
            line["lc_g_fr"] = linecode["g_fr"]
        end
    end
end


"""
    add_linecode_sequences!(data_eng)

Add linecode sequences to the `data_eng` dictionary.

# Arguments
- `data_eng`: A dictionary containing network data.

This function reads linecode data from a file and adds the linecode sequences to the `data_eng` dictionary.

"""
function add_linecode_sequences!(data_eng)
    #test opening linecode file
    linecode_file = joinpath("networks", "network_$(ntw)", "Feeder_$(fdr)", "LineCode.txt")
    linecode_data = open(linecode_file) do file
        readlines(file)
    end

    for (LineCodeName, lineCodeDic) in data_eng["linecode"]
        # Find the linecode data for the given LineCodeName
        linecode_entry = findfirst(contains.(lowercase.(linecode_data), lowercase(LineCodeName)))

        
        linecode_properties = split(linecode_data[linecode_entry], " ")
        linecode_name = linecode_values = split.(linecode_properties[2], ".")[2]
        linecode_values = split.(linecode_properties[3:end], "=")
        linecode_dict = Dict{String, String}()
        for linecode in linecode_values
            key = linecode[1]
            value = linecode[2]
            linecode_dict[key] = value
        end

        lineCodeDic["R0"] = linecode_dict["R0"]
        lineCodeDic["X0"] = linecode_dict["X0"]
        lineCodeDic["C0"] = linecode_dict["C0"]
        lineCodeDic["R1"] = linecode_dict["R1"]
        lineCodeDic["X1"] = linecode_dict["X1"]
        lineCodeDic["C1"] = linecode_dict["C1"]
        lineCodeDic["nphase"] = linecode_dict["nphases"]
        lineCodeDic["Units"] = linecode_dict["Units"]
    end
end



"""
    add_feeders_cktbks!(eng)

Add feeders and circuit breakers to the given `eng` dictionary.

# Arguments
- `eng`: A dictionary representing the engineering data.

# Description
This function iterates over the lines in the `eng` dictionary and adds the lines that match the pattern `feederd+` to the `Feeders` dictionary, and the lines that match the pattern `cktbkd+` to the `Cktbk` dictionary. Finally, it adds the `Cktbk` and `Feeders` dictionaries to the `eng` dictionary.

"""


function add_feeders_cktbks!(eng)
    Feeders = Dict{}()
    Cktbk = Dict{}()
    for (l, line) in eng["line"]
        if occursin(r"feeder\d+", l)
            Feeders[l] = line;
        end
    end
    for (cb, cktbk) in eng["line"]
        if occursin(r"cktbk\d+", cb)
            Cktbk[cb] = cktbk;
        end
    end
    eng["Cktbk"] = Cktbk;
    eng["Feeders"] = Feeders;
end


"""
    plot_network(eng; plot_bus=true, overlayNetwork=false, save_fig=false, plotname="plot.png"::String, save_json=false, json_filename="Networkfile.json"::String)

This function plots the network based on the given network data.

## Arguments
- `eng`: A dictionary containing the network data.
- `plot_bus`: A boolean indicating whether to plot the buses. Default is `true`.
- `overlayNetwork`: A boolean indicating whether to overlay the network on an existing plot. Default is `false`.
- `save_fig`: A boolean indicating whether to save the plot as an image file. Default is `false`.
- `plotname`: A string specifying the name of the image file to be saved. Default is "plot.png".
- `save_json`: A boolean indicating whether to save the network data as a JSON file. Default is `false`.
- `json_filename`: A string specifying the name of the JSON file to be saved. Default is "Networkfile.json".

## Returns
- `p`: The plot object representing the network plot.


"""
function plot_network(eng; plot_bus=true, overlayNetwork=false, save_fig=false, plotname="plot.png"::String, save_json=false, json_filename="Networkfile.json"::String)
    transformer_id = string.(keys(eng["transformer"])...)
    # extract feederno from eng["name"]
    feederno = parse(Int, match(r"\d+", eng["name"]).match)
    lats = []; lons = []; bus_labels = []; markershapes = []; markersizes = []; markercolors= []
    if plot_bus 
        
        for (bus_key, bus_data) in pairs(eng["bus"])
            if haskey(bus_data, "lat")
                push!(lats, bus_data["lat"])
                push!(lons, bus_data["lon"])
                push!(bus_labels, bus_key)
            else
                @warn("No coords for bus $bus_key")
            end
        end
        bus_hover_text = ["Bus: $bus_label" for bus_label in bus_labels] 
        

        three_phase_buses = [load["bus"] for (lo, load) in pairs(eng["load"]) if length(load["connections"]) > 2]
        phase_a_buses = [load["bus"] for (lo, load) in pairs(eng["load"]) if (length(load["connections"]) == 2 && load["connections"][1] == 1)]
        phase_b_buses = [load["bus"] for (lo, load) in pairs(eng["load"]) if (length(load["connections"]) == 2 && load["connections"][1] == 2)]
        phase_c_buses = [load["bus"] for (lo, load) in pairs(eng["load"]) if (length(load["connections"]) == 2 && load["connections"][1] == 3)]


        for (bus_key, bus_data) in pairs(eng["bus"])
            if bus_key in three_phase_buses
                push!(markershapes, :diamond)
                push!(markersizes, 5)
                push!(markercolors, "purple")
            elseif bus_key in phase_a_buses
                push!(markershapes, :dtriangle)
                push!(markersizes, 5)
                push!(markercolors, "red")
            elseif bus_key in phase_b_buses
                push!(markershapes, :dtriangle)
                push!(markersizes, 5)
                push!(markercolors, "green")
            elseif bus_key in phase_c_buses
                push!(markershapes, :dtriangle)
                push!(markersizes, 5)
                push!(markercolors, "blue")
            elseif bus_key == eng["Cktbk"]["cktbk$feederno"]["t_bus"]
                push!(markershapes, :octagon)
                push!(markersizes, 10)
                push!(markercolors, "yellow")
            else
                push!(markershapes, :circle)
                push!(markersizes, 1)
                push!(markercolors, "black")
            end
        end
        
    end 

    lines_x = []; lines_y = [] ; line_labels = [] ;line_colors = [] ; line_styles = []

    for (line_key, line_data) in pairs(eng["line"])
        f_bus_key = line_data["f_bus"]
        t_bus_key = line_data["t_bus"]
        if haskey(eng["bus"], f_bus_key) && haskey(eng["bus"], t_bus_key)
            if haskey(eng["bus"][f_bus_key], "lat") && haskey(eng["bus"][t_bus_key], "lat")
                f_bus_lat = eng["bus"][f_bus_key]["lat"]
                f_bus_lon = eng["bus"][f_bus_key]["lon"]
                t_bus_lat = eng["bus"][t_bus_key]["lat"]
                t_bus_lon = eng["bus"][t_bus_key]["lon"]
                push!(lines_x, [f_bus_lon, t_bus_lon])
                push!(lines_y, [f_bus_lat, t_bus_lat])
                push!(line_labels, "Line: $line_key: $f_bus_key to $t_bus_key")
                # Check if the line belongs to feeder1 and set its color accordingly
                if haskey(eng["line"][line_key], "Feeder")
                    push!(line_colors, "blue")
                else
                    push!(line_colors, "red")
                end


                if haskey(eng["line"][line_key], "Feeder")
                    push!(line_colors, "black")
                else
                    push!(line_colors, "red")
                end

                if occursin(r"feeder\d+",line_key)
                    push!(line_styles, :dot);
                elseif occursin(r"cktbk\d+",line_key)
                    push!(line_styles, :dash);
                else
                    push!(line_styles, :solid)
                end 

            end
        end
    end

    # Create a line plot for lines between buses with colors
    if overlayNetwork
        p = plot!(lines_x, lines_y, color=line_colors, hover=line_labels,
        linealpha=0.5, linewidth=1, 
        linestyle = line_styles, 
        legend=false)
    else
        p = plot(lines_x, lines_y, color=line_colors, hover=line_labels,
        linealpha=0.5, linewidth=1, size=(1000, 800),
        linestyle = line_styles,
        legend=false)
        
    end

    if plot_bus
    scatter!(p, lons, lats, color=markercolors, markersize=markersizes, markershape= markershapes, hover=bus_hover_text,
                label=markershapes, legend=false,
                axis = true, framestyle= :box, size=(1000, 800), title="Spanish Network, Transformer: $(transformer_id...) Feeder $feederno",
                )
    end
    if save_fig
        Plots.gr()
        savefig(p, plotname)
    end
    if save_json
        print_file(json_filename, eng)
    end
    return p
end



"""
    clean_irrelevant_data!(eng_feeder, transformer_id, feederno, flattened_member_lines)

This function cleans the irrelevant data from the `eng_feeder` dictionary based on the provided `transformer_id`, `feederno`, and `flattened_member_lines`.

## Arguments
- `eng_feeder`: A dictionary representing the feeder data.
- `transformer_id`: The ID of the transformer to keep in the `eng_feeder` dictionary.
- `feederno`: The feeder number to be used in the `eng_feeder` dictionary.
- `flattened_member_lines`: An array of line IDs to keep in the `eng_feeder` dictionary.



"""

function clean_irrelevant_data!(eng_feeder, transformer_id, feederno, flattened_member_lines)
    

    for (tr, transformer) in eng_feeder
        if tr != transformer_id
            delete!(eng_feeder["transformer"], tr)
        end 
    end

    for l in flattened_member_lines
        eng_feeder["line"][l]["Feeder"] = "feeder$feederno" 
    end


    # deleted unwanted lines
    for (l,lines) in eng_feeder["line"]

        if !in(l, flattened_member_lines)
            delete!(eng_feeder["line"], l)
        end
    end

    # delete unwanted buses
    keep_buses = []
    for (l,lines) in eng_feeder["line"]

        related_buses = [lines["f_bus"], lines["t_bus"]]
        for b in related_buses
            if !in(b, keep_buses)
                push!(keep_buses, b)
            end
        end

    end

    for (b,buses) in eng_feeder["bus"]
        if !in(b, keep_buses)
            delete!(eng_feeder["bus"], b)
        end
    end

    # delete unwanted loads
    for (lo, load) in eng_feeder["load"]
        if !in(load["bus"], keep_buses)
            delete!(eng_feeder["load"], lo)
        end
    end


    # delete unwanted transformers
    for (tr, transformer) in eng_feeder["transformer"]
        if tr != transformer_id[1]
            delete!(eng_feeder["transformer"], tr)
        end 
    end

    # delete unwanted linecodes
    keep_codes = []
    for (l,lines) in eng_feeder["line"]

        related_codes = [lines["linecode"]]
        for c in related_codes
            if !in(c, keep_codes)
                push!(keep_codes, c)
            end
        end

    end

    for (c,codes) in eng_feeder["linecode"]
        if !in(c, keep_codes)
            delete!(eng_feeder["linecode"], c)
        end
    end
    # delete unwanted Feeders
    for (f, feeder) in eng_feeder["Feeders"]
        if f != "feeder$feederno"
            delete!(eng_feeder["Feeders"], f)
        end
    end

    # delete unwanted CBs
    for (c, cktbk) in eng_feeder["Cktbk"]
        if c != "cktbk$feederno"
            delete!(eng_feeder["Cktbk"], c)
        end
    end


    # delete unwanted shunts

    for (sh, shunts) in eng_feeder["shunt"]
        if sh != "line_loop.grnd$feederno"
            delete!(eng_feeder["shunt"], sh)
        end
    end

    # delete unwanted loadshapes

    if haskey(eng_feeder, "time_series")
        keep_shapes = []
        
        for (lo, load) in eng_feeder["load"]
            related_shapes = [load["time_series"]["pd_nom"]]
            for s in related_shapes
                if !in(s, keep_shapes)
                    push!(keep_shapes, s)
                end
            end
        end
        
        for (sh, shape) in eng_feeder["time_series"]
            if !in(sh, keep_shapes)
                delete!(eng_feeder["time_series"], sh)
            end
        end
    end
        
    eng_feeder["name"] = "feeder$feederno"

end




"""
isolate_feeder(eng, feederno)

Given a network data structure `eng` and a feeder number `feederno`, this function creates a deep copy of the `eng` structure and isolates the specified feeder.

# Arguments
- `eng`: A network data structure.
- `feederno`: The feeder number to isolate.

# Returns
- `eng_feeder`: A deep copy of the `eng` structure with the specified feeder isolated.

# Description
The function isolates the specified feeder by performing the following steps:
1. Creates a deep copy of the `eng` structure.
2. Sets the voltage source bus to the TrLVBus of the specified feeder.
3. Retrieves the transformer ID(s) associated with the TrLVBus.
4. Retrieves the feeder head bus and the circuit breaker buses of the specified feeder.
5. Constructs a list of member lines by traversing the network starting from the feeder head bus.
6. Filters out any empty member lines and flattens the list.
7. Appends the transformer ID(s), feeder number, and circuit breaker number to the list of member lines.
8. Removes irrelevant data from the `eng_feeder` structure based on the transformer ID(s), feeder number, and member lines.

"""
function isolate_feeder(eng, feederno)
    eng_feeder  = deepcopy(eng)


    TrLVBus = eng_feeder["Feeders"]["feeder$feederno"]["f_bus"]
    eng_feeder["voltage_source"]["source"]["bus"] = TrLVBus 


    transformer_id = [x for (x, xfmr) in pairs(eng_feeder["transformer"]) if xfmr["bus"][2] == TrLVBus]
    display("Network: $(transformer_id...), Feeder: $feederno")
    Feeder_head = eng_feeder["Feeders"]["feeder$feederno"]["t_bus"]
    Feeder_head = eng_feeder["Cktbk"]["cktbk$feederno"]["f_bus"]
    cb_2 = eng_feeder["Cktbk"]["cktbk$feederno"]["t_bus"]

    ################
    moving_LVBus = [Feeder_head]
    member_lines = []
    while length(moving_LVBus) > 0
        temp_lines = []
        while length(moving_LVBus) > 0
            member_line = [x for (x, line) in pairs(eng_feeder["line"]) if line["f_bus"] == moving_LVBus[end]]
            push!(member_lines, member_line)
            push!(temp_lines, member_line...)
            pop!(moving_LVBus)
        end
        moving_LVBus = [string(eng_feeder["line"][l]["t_bus"]) for l in temp_lines]
    end
    ##################


    flattened_member_lines = filter(x -> length(x) > 0, vcat(member_lines...))
    mv_id= string.([x for (x, line) in eng_feeder["line"] if line["t_bus"] == eng_feeder["transformer"][string.(transformer_id...)]["bus"][1]]...)
    push!(flattened_member_lines, mv_id) 
    push!(flattened_member_lines, "feeder$feederno")
    push!(flattened_member_lines, "cktbk$feederno")
    
 
    clean_irrelevant_data!(eng_feeder, transformer_id, feederno, flattened_member_lines)    
    return eng_feeder
end


"""
    attach_Yg_ground!(eng)

This function attaches the Yg ground to the specified transformer in the given `eng` structure.

# Arguments
- `eng`: The structure containing the network data.

# Description
- The function sets the base frequency to 50 and the default sbase to 22000.
- For each transformer in the `eng` structure, it finds the grounded bus.
- It then calculates the grnd_react value based on the line data.
- The function updates the `xg` and `rg` values of the grounded bus using the `grnd_react` value.
- Finally, it deletes the line with the `grnd_react` value from the `eng` structure.

"""
function attach_Yg_ground!(eng)
    eng["settings"]["base_frequency"] = 50  # because it is was 60  
    eng["settings"]["sbase_default"] = 22000 # because it was 100000
    for (tr, transformer) in eng["transformer"]
        grounded_bus = transformer["bus"][2]
        grnd_react = string.([x for (x, line) in eng["line"] if (line["f_bus"] == grounded_bus && haskey(line, "xs"))]...)
        eng["bus"][grounded_bus]["xg"] = eng["line"][grnd_react]["xs"]
        eng["bus"][grounded_bus]["rg"] = eng["line"][grnd_react]["rs"]
        delete!(eng["line"], grnd_react)
    end
end

"""
    add_vcomplex_vm_va!(res)

This function calculates the complex voltage, magnitude, and angle for each bus in the given network data.

# Arguments
- `res`: A dictionary containing the network data.

"""
function add_vcomplex_vm_va!(res)


    if haskey(res["solution"], "nw")           

        for (nw, network) in res["solution"]["nw"]
            for (b, bus) in network["bus"]
                bus["vcomplex"] = bus["vr"] + bus["vi"]*im
                bus["vm"] = abs.(bus["vcomplex"])
                bus["va"] = angle.(bus["vcomplex"])*180/π
            end
        end

    else

     for (b, bus) in res["solution"]["bus"]
                bus["vcomplex"] = bus["vr"] + bus["vi"]*im
                bus["vm"] = abs.(bus["vcomplex"])
                bus["va"] = angle.(bus["vcomplex"])*180/π
     end
    
    end 

end



function plot_bus_ts_vm(res, bus, Hours)
    
    
    vm_a = []
    vm_b = []
    vm_c = []
    vm_n = []

    for (_, network) in res["solution"]["nw"]
        push!(vm_a, network["bus"]["$bus"]["vm"][1])
        push!(vm_b, network["bus"]["$bus"]["vm"][2])
        push!(vm_c, network["bus"]["$bus"]["vm"][3])
        push!(vm_n, network["bus"]["$bus"]["vm"][4])
    end

    timeseries_df_vm = DataFrame(network = Int[], vm_a = Float64[], vm_b = Float64[], vm_c = Float64[], vm_n = Float64[])

    for i in 1:Hours
        push!(timeseries_df_vm, (i, vm_a[i], vm_b[i], vm_c[i], vm_n[i]))
    end


    plotlyjs()
    plot_abc = plot(x=timeseries_df_vm.network, timeseries_df_vm.vm_a, label="A", xlabel="Time (Hours)", ylabel="Voltage Magnitude (pu)");
    plot!(plot_abc, x=timeseries_df_vm.network, timeseries_df_vm.vm_b, label="B");
    plot!(plot_abc, x=timeseries_df_vm.network, timeseries_df_vm.vm_c, label="C");
    plot_n = plot(x=timeseries_df_vm.network, timeseries_df_vm.vm_n, label="N", xlabel="Time (Hours)", ylabel="Voltage Magnitude (pu)");

    # create a layout 

    vm_plot = plot(plot_abc, plot_n, layout=(2, 1), legend=:topleft, link=:x,
        linestyle=[:solid :solid :solid :solid],
        linecolor=[:red :blue :green :black],
        title="Bus $bus Voltage Magnitude",
        fontfamily="Times",    
        xticks = 2:2:Hours,
        #yticks = 0:0.0001:1.1,
        size=(1000,700),
        )


    return vm_plot
end



function plot_bus_ts_vm(p, res, bus, Hours)
    
    
    vm_a = []
    vm_b = []
    vm_c = []
    vm_n = []

    for (_, network) in res["solution"]["nw"]
        push!(vm_a, network["bus"]["$bus"]["vm"][1])
        push!(vm_b, network["bus"]["$bus"]["vm"][2])
        push!(vm_c, network["bus"]["$bus"]["vm"][3])
        push!(vm_n, network["bus"]["$bus"]["vm"][4])
    end

    timeseries_df_vm = DataFrame(network = Int[], vm_a = Float64[], vm_b = Float64[], vm_c = Float64[], vm_n = Float64[])

    for i in 1:Hours
        push!(timeseries_df_vm, (i, vm_a[i], vm_b[i], vm_c[i], vm_n[i]))
    end


    plotlyjs()
    plot_abc = plot(x=timeseries_df_vm.network, timeseries_df_vm.vm_a, label="A", xlabel="Time (Hours)", ylabel="Voltage Magnitude (pu)");
    plot!(plot_abc, x=timeseries_df_vm.network, timeseries_df_vm.vm_b, label="B");
    plot!(plot_abc, x=timeseries_df_vm.network, timeseries_df_vm.vm_c, label="C");
    plot_n = plot(x=timeseries_df_vm.network, timeseries_df_vm.vm_n, label="N", xlabel="Time (Hours)", ylabel="Voltage Magnitude (pu)");

    # create a layout 

    vm_plot = plot!(p, plot_abc, plot_n, layout=(2, 1), legend=:topleft, link=:x,
        linestyle=[:solid :solid :solid :solid],
        linecolor=[:red :blue :green :black],
        title="Bus $bus Voltage Magnitude",
        fontfamily="Times",    
        xticks = 2:2:Hours,
        #yticks = 0:0.0001:1.1,
        size=(1000,700),
        )


    return vm_plot
end









function plot_bus_ts_va(res, bus, Hours)

    va_a = []
    va_b = []
    va_c = []
    va_n = []

    for (_, network) in res["solution"]["nw"]
        push!(va_a, network["bus"]["$bus"]["va"][1])
        push!(va_b, network["bus"]["$bus"]["va"][2])
        push!(va_c, network["bus"]["$bus"]["va"][3])
        push!(va_n, network["bus"]["$bus"]["va"][4])
    end

    timeseries_df_va = DataFrame(network = Int[], va_a = Float64[], va_b = Float64[], va_c = Float64[], va_n = Float64[])
    for i in 1:Hours
        push!(timeseries_df_va, (i, va_a[i], va_b[i], va_c[i], va_n[i]))
    end


    plot_va_a = plot(x=timeseries_df_va.network, timeseries_df_va.va_a, label="A", xlabel="Time (Hours)", ylabel="Voltage Angle (degrees)");
    plot_va_b= plot(x=timeseries_df_va.network, timeseries_df_va.va_b, label="B", xlabel="Time (Hours)", ylabel="Voltage Angle (degrees)");
    plot_va_c = plot(x=timeseries_df_va.network, timeseries_df_va.va_c, label="C", xlabel="Time (Hours)", ylabel="Voltage Angle (degrees)");
    plot_va_n = plot(x=timeseries_df_va.network, timeseries_df_va.va_n, label="N", xlabel="Time (Hours)", ylabel="Voltage Angle (degrees)");

    # create a layout

    va_plot = plot(plot_va_a, plot_va_b, plot_va_c, plot_va_n, layout=(4, 1), legend=:topleft, link=:x,
        linestyle=[:solid :solid :solid :solid],
        linecolor=[:red :blue :green :black],
        title="Bus $bus Voltage Angle",
        fontfamily="Times",    
        xticks = 2:2:Hours,
        yticks = -180:0.0001:180,
        size=(1000,700),
        )

    return va_plot 

end










function plot_load_ts_cid(res, lo, Hours)

    cid = []

    for (_, network) in res["solution"]["nw"]
        push!(cid, -network["load"]["$lo"]["cid"][1]*res["solution"]["nw"]["1"]["settings"]["sbase_default"])
    end

    timeseries_df_cid = DataFrame(network = Int[], cid = Float64[])
    for i in 1:Hours
        push!(timeseries_df_cid, (i, cid[i]))
    end
    plot_cid = plot(x=timeseries_df_cid.network, timeseries_df_cid.cid,
                     label="CID",
                     xlabel="Time (Hours)",
                     ylabel="Power (W)",
                     fontfamily="Times",    
                     xticks = 2:2:Hours,
                     size=(1000,350),
                     framestyle = :box,
                     
                     
                     
                     );
    return plot_cid    
end






function plot_load_ts_cid(p, res, lo, Hours)

    cid = []

    for (_, network) in res["solution"]["nw"]
        push!(cid, network["load"]["$lo"]["cid"][1])
    end

    timeseries_df_cid = DataFrame(network = Int[], cid = Float64[])
    for i in 1:Hours
        push!(timeseries_df_cid, (i, cid[i]))
    end
    plot_cid = plot!(p, x=timeseries_df_cid.network, timeseries_df_cid.cid,
                     label="CID",
                     xlabel="Time (Hours)",
                     ylabel="Power (pu)",
                     fontfamily="Times",    
                     xticks = 2:2:Hours,
                     size=(1000,350),
                     framestyle = :box,
                     
                     
                     
                     );
    return plot_cid    
end






"""
    merge_bus_results!(res_mn, mn_Feeder_math)

Merge the bus results from `res_mn` into `mn_Feeder_math`.

# Arguments
- `res_mn`: A dictionary containing the bus results.
- `mn_Feeder_math`: A dictionary representing the network data.

# Description
This function merges the bus results from `res_mn` into the `mn_Feeder_math` dictionary. It updates the `"vi"`, `"va"`, and `"vcomplex"` values for each bus in the network.

"""
function merge_bus_results!(res_mn, mn_Feeder_math)
    for (nw, _) in mn_Feeder_math["nw"]
        for (bus, bus_data) in mn_Feeder_math["nw"][nw]["bus"]
            bus_data["vi"] = res_mn["solution"]["nw"][nw]["bus"][bus]["vi"]
            bus_data["va"] = res_mn["solution"]["nw"][nw]["bus"][bus]["va"]
            bus_data["vcomplex"] = res_mn["solution"]["nw"][nw]["bus"][bus]["vcomplex"]
        end
    end
end




"""
    attach_load_details!(res_mn, mn_Feeder_math)

This function attaches load details from the `res_mn` solution to the `mn_Feeder_math` data structure.

**Arguments**
- `res_mn`: A dictionary containing the solution results.
- `mn_Feeder_math`: A data structure representing the feeder network.

**Returns**
- Nothing

"""


function attach_load_details!(res_mn, mn_Feeder_math)

    for (nw, _) in mn_Feeder_math["nw"]
        for (lo, load) in mn_Feeder_math["nw"][nw]["load"]
            load["qd_bus"] = res_mn["solution"]["nw"][nw]["load"][lo]["qd_bus"]
            load["pd_bus"] = res_mn["solution"]["nw"][nw]["load"][lo]["pd_bus"]
            load["qd"] = res_mn["solution"]["nw"][nw]["load"][lo]["qd"]
            load["pd"] = res_mn["solution"]["nw"][nw]["load"][lo]["pd"]
            load["cid_bus"] = res_mn["solution"]["nw"][nw]["load"][lo]["cid_bus"]
            load["cid"] = res_mn["solution"]["nw"][nw]["load"][lo]["cid"]
            load["crd_bus"] = res_mn["solution"]["nw"][nw]["load"][lo]["crd_bus"]
            load["crd"] = res_mn["solution"]["nw"][nw]["load"][lo]["crd"]
        end
    end
    
    
    for (nw, _) in mn_Feeder_math["nw"]
        for (lo, load) in mn_Feeder_math["nw"][nw]["load"]
            load["vi"] = res_mn["solution"]["nw"][nw]["bus"][string.(load["load_bus"])]["vi"]
            load["vr"] = res_mn["solution"]["nw"][nw]["bus"][string.(load["load_bus"])]["vr"]
            load["va"] = res_mn["solution"]["nw"][nw]["bus"][string.(load["load_bus"])]["va"]
            load["vm"] = res_mn["solution"]["nw"][nw]["bus"][string.(load["load_bus"])]["vm"]
            load["vcomplex"] = res_mn["solution"]["nw"][nw]["bus"][string.(load["load_bus"])]["vcomplex"]
        end
    end
    
    
end


function write_dss_file(spanish_data::Dict, feeder_no::Int) 
    Feeder = isolate_feeder(spanish_data, feeder_no)

    transformer_name = first(Feeder["transformer"]).second["name"]
    SecondaryV = first(Feeder["transformer"]).second["vm_nom"][2]
    PrimaryBus = first(Feeder["transformer"]).second["bus"][1]
    SecondaryBus = first(Feeder["transformer"]).second["bus"][2]

    line_bf = string.([x for(x, line) in Feeder["line"] if line["f_bus"] == PrimaryBus || line["t_bus"] == PrimaryBus]...)
    delete!(Feeder["line"], line_bf)
    dss_file = open("C:/Users/mnumair/OneDrive - KU Leuven/PhD Agenda/2024/24-11 November/TSK-478 mappingEULV/Mapping Spanish/spanish_feeders/Feeder_$(feeder_no)_Ntw_$(transformer_name).dss", "w")

    # WRITING THE DSS FILE
    println(dss_file, "clear")
    println(dss_file, "new circuit.Feeder_$(feeder_no)_Ntw_$(transformer_name)")
    println(dss_file, "edit Vsource.Source bus1=$SecondaryBus basekv=$SecondaryV pu=1 phases=3 ISC3=9000  ISC1=5000")
    println(dss_file, "! edit Vsource.Source bus1=source basekv=22 pu=1.0001 phases=3 ISC3=9000  ISC1=5000 ! only uncomment if you want to include the transfomer")
    
    println(dss_file, "!-- Linecodes ------------------------------------------------------------------------------------------------------")
    println(dss_file, "redirect linecode.txt")
    
    println(dss_file, "!-- Loadshapes ------------------------------------------------------------------------------------------------------")
    println(dss_file, "redirect Loadshape.txt")
    
    println(dss_file, "! -TRFRMS-----------------------------------------------------------------------------------------------------------------")
    
    for (tr, transformer) in Feeder["transformer"]
        #buses = [transformer["bus"][i] * "." * join(transformer["connections"][i], ".") for i in 1:2]
        # buses = [transformer["bus"][1],
        #          transformer["bus"][2] * "." * join([1 2 3 4], ".")]
        
        println(dss_file, "! The transformer is commented out for an LV feeder analysis. You can uncomment but uncomment the Vsource line above")
        println(dss_file, "! New line.mv Bus1=sourcebus.1.2.3 Bus2=$(transformer["bus"][1]).1.2.3 phases=3 Linecode=101 Length=5 Units=m")
        
        
        buses = [transformer["bus"][1],
                 transformer["bus"][2] * "." * join([1 2 3 4], ".")]
    
    
    
        conns = [string(transformer["configuration"][i]) for i in 1:2]
        kVs = [transformer["vm_nom"][i] for i in 1:2]
        kVAs = [transformer["sm_nom"][i] for i in 1:2]
        XHL = transformer["xsc"][1]
    
        println(dss_file, "! New Transformer." * transformer["name"] *
                " windings=2"*
                " Buses=[" * join(buses, " ") * "]" *
                " Conns=[" * join(conns, " ") * "]" *
                " kVs=[" * join(string.(kVs), " ") * "]" *
                " kVAs=[" * join(string.(kVAs), " ") * "]" *
                " XHL=" * string(XHL*100) *
                " phases=3"*
                " sub=y ")
            
            
            println(dss_file, "New Reactor.Yngrnd phases=1 bus1=$(transformer["bus"][2]).4 bus2=$(transformer["bus"][2]).0 R=5 X=0.01 ")
        end
    
    
    println(dss_file, "! -LINES-----------------------------------------------------------------------------------------------------------------")
    
    for (l, line) in Feeder["line"]
        f_phases = join(line["f_connections"], ".")
        t_phases = join(line["t_connections"], ".")
        println(dss_file, "New Line." * l * 
                " Bus1=" * string(line["f_bus"]) * "." * f_phases * 
                " Bus2=" * string(line["t_bus"]) * "." * t_phases * 
                " phases=$(length(line["t_connections"]))" * 
                " linecode=" * string(line["linecode"]) * 
                " Length=" * string(line["length"]) * 
                " Units=m")
    end
    
    println(dss_file, "! --LOADS-----------------------------------------------------------------------------------------------------------")
    for (lo, load) in Feeder["load"]
        v_mult = length(load["connections"]) > 2 ? sqrt(3) : 1
        bus = load["bus"] * "." * join(load["connections"], ".")
        kV = load["vm_nom"]*v_mult
        kW = sum(load["pd_nom"])  # Sum the power for all phases
        PF = 0.95  # Assuming Power Factor is constant as 0.95
        daily = load["time_series"]["pd_nom"]
        phases = length(load["connections"])-1
    
        println(dss_file, "New " * load["source_id"] * 
                " Phases=" * string(phases) * 
                " Bus1=" * bus * 
                " kV=" * string(kV) * 
                " kW=" * string(kW) * 
                " PF=" * string(PF) * 
                " daily=" * daily)
    end
    
    
    println(dss_file, "! --PF -----------------------------------------------------------------------------------------------------------
    set mode=Daily number=481 stepsize=60m
    Calcvoltagebases
    solve
    BusCoords BusCoords_ETRS.csv
    Plot circuit Power Max=1000 dots=n labels=n subs=y C1=black 1ph=3
    ")

    close(dss_file)

end 


# Coord Transform

transSpanishTo4326 = Proj.Transformation("EPSG:3042", "EPSG:4326", always_xy=true) # ETRS89 / UTM zone 30N to WGS84
#transSpanishTo4326 = Proj.Transformation("EPSG:25830", "EPSG:4326", always_xy=true) # ETRS89 / UTM zone 30N to WGS84


function extract_spanish_network_and_feeder(file_path::String; pattern = r"(?i)\\Feeder_(\d+)_Ntw_(\d+).dss")
    
    # Match the pattern against the file path
    matchedPattern = match(pattern, file_path)
    
    if matchedPattern !== nothing
        # Extract the network and feeder numbers
        Fdr = parse(Int, matchedPattern.captures[1])
        Ntw = parse(Int, matchedPattern.captures[2])
        return Ntw, Fdr
    else
        error("Pattern not found in the given file path")
    end
end