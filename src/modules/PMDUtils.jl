"""
    PMDUtils

Internal sub-module providing utility functions for PowerModelsDistribution (PMD) workflows.

This module re-exports functions from the main Pliers module for:
- Processing power flow solutions (voltage fluffing, dictification)
- Impedance calculations (Kron reduction, sequence components)
- Network data manipulation
- Result analysis and transformation
- Engineering and mathematical model exploration

See the main Pliers module for function documentation.
"""
module PMDUtils

using ..Pliers



#=
░██████████ ░███    ░██   ░██████     ░█████████  ░██████████ ░█████████    ░██████   ░█████████  ░██████████
░██         ░████   ░██  ░██   ░██    ░██     ░██ ░██         ░██     ░██  ░██   ░██  ░██     ░██     ░██    
░██         ░██░██  ░██ ░██           ░██     ░██ ░██         ░██     ░██ ░██     ░██ ░██     ░██     ░██    
░█████████  ░██ ░██ ░██ ░██  █████    ░█████████  ░█████████  ░█████████  ░██     ░██ ░█████████      ░██    
░██         ░██  ░██░██ ░██     ██    ░██   ░██   ░██         ░██         ░██     ░██ ░██   ░██       ░██    
░██         ░██   ░████  ░██  ░███    ░██    ░██  ░██         ░██          ░██   ░██  ░██    ░██      ░██    
░██████████ ░██    ░███   ░█████░█    ░██     ░██ ░██████████ ░██           ░██████   ░██     ░██     ░██    



=#



## eng_report
"""
    eng_report(eng::Dict{String, Any})

Generate a report for the electrical network described by the dictionary `eng`.

# Arguments
- `eng::Dict{String, Any}`: A dictionary containing various components of the electrical network.

# Description
This function extracts various components from the `eng` dictionary, such as buses, lines, linecodes, loads, voltage sources, time series data, conductor IDs, name, settings, files, and data model. It then prints a formatted report summarizing the contents of the network, including the number of each component present.

# Example
```julia
using PowerModelsDistribution
using Pliers 
eng= PowerModelsDistribution.parse_file("example.dss")
eng_report(eng)
```
"""
function eng_report(eng::Dict{String,Any}; detailed=false)
    print(UNDERLINE(BLUE_FG("Report for the ", BOLD("$(get(eng, "data_model", nothing))"),
        " model of the network  ", BOLD("$(get(eng, "name", nothing))"), " prased from ", BOLD("$(split(get(eng,"files",[nothing])[1],"/")[end])"), " \n")))

    print(WHITE_FG("This network has:\n"))
    print("             $(BOLD("$(length(get(eng, "voltage_source", [])))")) voltage sources, \n")
    print("             $(BOLD("$(length(get(eng, "bus", [])))")) buses, \n")
    print("             $(BOLD("$(length(get(eng, "conductor_ids", [])))")) terminals per bus, \n")
    print("             $(BOLD("$(length(get(eng, "line", [])))")) lines, \n")
    print("             $(BOLD("$(length(get(eng, "linecode", [])))")) linecodes, \n")
    print("             $(BOLD("$(length(get(eng, "load", [])))")) loads, \n")
    print("             $(BOLD("$(length(get(eng, "transformer", [])))")) transformers, \n")
    print("             $(BOLD("$(length(get(eng, "time_series", [])))")) time series data points.\n")

    _eng_summary_table(eng)

    if detailed == true
        buses_table(eng)
        lines_table(eng)
        loads_table(eng)
        if haskey(eng, "transformer") && length(eng["transformer"]) > 0
            transformers_table(eng)
        end
    end

end

function _eng_summary_table(eng::Dict{String,Any})
    n_nodes = length(get(eng, "bus", Dict()))
    n_lines = length(get(eng, "line", Dict()))
    n_loads = length(get(eng, "load", Dict()))
    n_tfs = length(get(eng, "transformer", Dict()))

    tot_len = 0.0
    tot_r = 0.0
    tot_x = 0.0
    max_amps = 0.0

    lines = get(eng, "line", Dict())
    linecodes = get(eng, "linecode", Dict())

    for (_, line) in lines
        l_len = get(line, "length", 0.0)
        tot_len += l_len

        rating = get(line, "normamps", get(line, "rate_a", 0.0))
        if isa(rating, AbstractArray) && !isempty(rating)
            rating = maximum(rating)
        elseif isa(rating, Number)
            # ok
        else
            rating = 0.0
        end
        if rating > max_amps
            max_amps = rating
        end

        r_line = 0.0
        x_line = 0.0

        # Resistance
        if haskey(line, "Rs (Ω)")
            val = line["Rs (Ω)"]
            r_line = isa(val, AbstractArray) ? val[1] : val
        elseif haskey(line, "rs") && !haskey(line, "linecode")
            val = line["rs"]
            r_val = isa(val, AbstractArray) ? val[1] : val
            r_line = r_val * l_len
        elseif haskey(line, "linecode")
            lc_id = line["linecode"]
            if haskey(linecodes, lc_id)
                lc = linecodes[lc_id]
                rs = get(lc, "rs", 0.0)
                r_val = isa(rs, AbstractArray) ? rs[1] : rs
                r_line = r_val * l_len
            end
        end
        tot_r += r_line

        # Reactance
        if haskey(line, "Xs (Ω)")
            val = line["Xs (Ω)"]
            x_line = isa(val, AbstractArray) ? val[1] : val
        elseif haskey(line, "xs") && !haskey(line, "linecode")
            val = line["xs"]
            x_val = isa(val, AbstractArray) ? val[1] : val
            x_line = x_val * l_len
        elseif haskey(line, "linecode")
            lc_id = line["linecode"]
            if haskey(linecodes, lc_id)
                lc = linecodes[lc_id]
                xs = get(lc, "xs", 0.0)
                x_val = isa(xs, AbstractArray) ? xs[1] : xs
                x_line = x_val * l_len
            end
        end
        tot_x += x_line
    end

    tot_mva = 0.0
    for (_, tf) in get(eng, "transformer", Dict())
        sm = get(tf, "sm_nom", 0.0)
        val = isa(sm, AbstractArray) && !isempty(sm) ? maximum(sm) : (isa(sm, Number) ? sm : 0.0)
        tot_mva += val
    end
    tot_mva_val = tot_mva / 1000.0

    lats = [b["lat"] for (k, b) in get(eng, "bus", Dict()) if haskey(b, "lat")]
    lons = [b["lon"] for (k, b) in get(eng, "bus", Dict()) if haskey(b, "lon")]
    area = 0.0
    gis_info = "No"
    if !isempty(lats)
        gis_info = "Yes"
        min_lat, max_lat = minimum(lats), maximum(lats)
        min_lon, max_lon = minimum(lons), maximum(lons)
        lat_dist = (max_lat - min_lat) * 111.0
        avg_lat = mean(lats)
        lon_dist = (max_lon - min_lon) * 111.0 * cos(deg2rad(avg_lat))
        area = lat_dist * lon_dist
    end

    rx = (tot_x != 0) ? tot_r / tot_x : 0.0
    z_tot = sqrt(tot_r^2 + tot_x^2)
    len_load = (n_loads > 0) ? tot_len / n_loads : 0.0
    z_cons = (n_loads > 0) ? z_tot / n_loads : 0.0
    n_linecodes = length(linecodes)

    data = [
        "Number of nodes" n_nodes;
        "Number of branches" n_lines;
        "Number of loads" n_loads;
        "Number Transformers" n_tfs;
        "Total resistance" round(tot_r, digits=4);
        "Total reactance" round(tot_x, digits=4);
        "Total MVA of transformers" round(tot_mva_val, digits=4);
        "Total length of lines (m)" round(tot_len, digits=4);
        #"Area (sq km)"              round(area, digits=4);
        "Number of line types used" n_linecodes;
        "Total impedance" round(z_tot, digits=4);
        "R/X ratio" round(rx, digits=4);
        "Length/load" round(len_load, digits=4);
        "Z/ consumer" round(z_cons, digits=4)
        #"GIS information"           gis_info;
        #"Max current capacity"      round(max_amps, digits=4);
    ]
    #display(data)
    println("\"Feeder name\", \"Number of nodes\", \"Number of branches\", \"Number of loads\", \"Number Transformers\", \"Total resistance\", \"Total reactance\", \"Total MVA of transformers\", \"Total length of lines (m)\", \"Number of line types used\", \"Total impedance\", \"R/X ratio\", \"Length/load\", \"Z/ consumer\"")
    println("$(eng["name"]),$(n_nodes), $(n_lines), $(n_loads), $(n_tfs), $(round(tot_r, digits=4)), $(round(tot_x, digits=4)), $(round(tot_mva_val, digits=4)), $(round(tot_len, digits=4)), $(n_linecodes), $(round(z_tot, digits=4)), $(round(rx, digits=4)), $(round(len_load, digits=4)), $(round(z_cons, digits=4))")
    #return data
    #pretty_table(data; header=["Attributes", "Values"], alignment=[:l, :r])
end


## buses_table

function _init_buses_df()
    return DataFrame(
        bus_id=String[],
        status=String[],
        terminals=Array{Int64,1}[],
        lat=Float64[],
        lon=Float64[],
        rg=Array{Float64,1}[],
        xg=Array{Float64,1}[],
        grounded=Array{Int64,1}[]
    )

end

function _push_buses_df!(buses_df, bus_id, bus)
    push!(buses_df, (
        bus_id,
        string(get(bus, "status", "")),
        get(bus, "terminals", Int[]),
        get(bus, "lat", 0.0),
        get(bus, "lon", 0.0),
        get(bus, "rg", Float64[]),
        get(bus, "xg", Float64[]),
        get(bus, "grounded", Int[])
    ))
end

"""
    buses_table(eng::Dict{String, Any})

Generate a table summarizing the buses in the electrical network described by the dictionary `eng`.

# Arguments
- `eng::Dict{String, Any}`: A dictionary containing various components of the electrical network.

# Description
This function extracts the buses from the `eng` dictionary and creates a DataFrame with the bus ID, status, terminals, resistance to ground (rg), reactance to ground (xg), and grounding status. It then prints a formatted table of the buses.

# Example
```julia
using PowerModelsDistribution
using Pliers 
eng= PowerModelsDistribution.parse_file("example.dss")
buses_table(eng)
```
"""
function buses_table(eng::Dict{String,Any})
    buses = haskey(eng, "bus") ? eng["bus"] : error("No buses found in the engineering data check the data model has \"bus\" key")
    buses_df = _init_buses_df()
    for (bus_id, bus) in buses
        _push_buses_df!(buses_df, bus_id, bus)
    end
    header("Buses Table ($(nrow(buses_df)) buses)")
    pretty_table(sort!(buses_df))
    extra_keys(buses, ["status", "terminals", "rg", "xg", "grounded", "lat", "lon"])
end

"""
    buses_table(eng::Dict{String,Any}, condition)

Generate and display a filtered table of buses from the given engineering data.

# Arguments
- `eng::Dict{String,Any}`: A dictionary containing engineering data, which must include a "bus" key with bus information.
- `condition`: A function that takes a bus dictionary as input and returns a boolean indicating whether the bus meets the filtering criteria.

# Description
This function extracts bus information from the provided engineering data dictionary, applies the given condition to filter the buses, and then displays the filtered buses in a formatted table. Each bus is augmented with its `bus_id` before filtering. The table includes columns for `bus_id`, `status`, `terminals`, `rg`, `xg`, and `grounded`.

# Example
```julia
using PowerModelsDistribution
using Pliers 
eng= Pliers.parse_file("example.dss")
buses_table(eng, bus -> bus["bus_id"] =="sourcebus")

````
or 

```julia
buses_table(eng, bus -> haskey(bus, "grounded") && bus["grounded"]==[4])
```
"""
function buses_table(eng::Dict{String,Any}, condition::Function)
    buses = haskey(eng, "bus") ? eng["bus"] : error("No buses found in the engineering data check the data model has \"bus\" key")
    # adding the bus_id to the bus dictionary so it can be filtered by name
    for (bus_id, bus) in buses
        bus["bus_id"] = bus_id
    end

    buses_df = _init_buses_df()
    matched_buses_idx = []

    for (bus_id, bus) in buses
        if condition(bus)
            _push_buses_df!(buses_df, bus_id, bus)
            push!(matched_buses_idx, bus_id)
        end
    end
    header("Filtered Buses Table ($(nrow(buses_df)) buses)")
    pretty_table(sort!(buses_df))
    extra_keys(buses, ["status", "terminals", "rg", "xg", "grounded", "lat", "lon"])
    return matched_buses_idx
end


## lines table

function _init_lines_df()
    return DataFrame(
        line_id=String[],
        source_id=String[],
        f_bus=String[],
        f_connections=Array{Int64,1}[],
        t_bus=String[],
        t_connections=Array{Int64,1}[],
        length=Float64[],
        linecode=String[],
        status=String[]
    )
end

function _push_lines_df!(lines_df, line_id, line)
    push!(lines_df, (
        line_id,
        get(line, "source_id", ""),
        get(line, "f_bus", ""),
        get(line, "f_connections", Int[]),
        get(line, "t_bus", ""),
        get(line, "t_connections", Int[]),
        get(line, "length", 0.0),
        get(line, "linecode", ""),
        string(get(line, "status", ""))
    ))
end


"""
    lines_table(eng::Dict{String, Any})

Generate a table summarizing the lines in the electrical network described by the dictionary `eng`.

# Arguments
- `eng::Dict{String, Any}`: A dictionary containing various components of the electrical network.

# Description
This function extracts the lines from the `eng` dictionary and creates a DataFrame with the line ID, status, from bus, to bus, length, resistance, reactance, and linecode. It then prints a formatted table of the lines.

# Example
```julia
using PowerModelsDistribution
using Pliers
eng= PowerModelsDistribution.parse_file("example.dss")
lines_table(eng)
```
"""
function lines_table(eng::Dict{String,Any})
    lines = haskey(eng, "line") ? eng["line"] : error("No lines found in the engineering data check the data model has \"line\" key")

    lines_df = _init_lines_df()

    for (line_id, line) in lines

        _push_lines_df!(lines_df, line_id, line)

    end
    header("Lines Table ($(nrow(lines_df)) lines)")
    pretty_table(sort!(lines_df, [:f_bus, :t_bus, :line_id]))

    extra_keys(lines, ["line_id", "source_id", "status", "f_bus", "f_connections", "t_bus", "t_connections", "length", "linecode"])

end

"""
    lines_table(eng::Dict{String,Any}, condition)

Generate and display a filtered table of lines from the given engineering data.

# Arguments
- `eng::Dict{String,Any}`: A dictionary containing engineering data, which must include a "line" key with line information.
- `condition`: A function that takes a line dictionary as input and returns a boolean indicating whether the line meets the filtering criteria.

# Description
This function extracts line information from the provided engineering data dictionary, applies the given condition to filter the lines, and then displays the filtered lines in a formatted table. Each line is augmented with its `line_id` before filtering. The table includes columns for `line_id`, `source_id`, `f_bus`, `f_connections`, `t_bus`, `t_connections`, `length`, `linecode`, and `status`.

# Example
```julia
using PowerModelsDistribution
using Pliers
eng= PowerModelsDistribution.parse_file("example.dss")
lines_table(eng, line -> line["length"] > 0.75)
```
"""
function lines_table(eng::Dict{String,Any}, condition::Function)
    lines = haskey(eng, "line") ? eng["line"] : error("No lines found in the engineering data check the data model has \"line\" key")
    # adding the line_id to the line dictionary so it can be filtered by name
    for (line_id, line) in lines
        line["line_id"] = line_id
    end

    lines_df = _init_lines_df()
    mathced_lines_idx = []
    for (line_id, line) in lines
        if condition(line)
            _push_lines_df!(lines_df, line_id, line)
            push!(mathced_lines_idx, line_id)
        end
    end
    header("Filtered Lines Table ($(nrow(lines_df)) lines)")
    pretty_table(sort!(lines_df, [:f_bus, :t_bus, :line_id]))
    extra_keys(lines, ["line_id", "source_id", "status", "f_bus", "f_connections", "t_bus", "t_connections", "length", "linecode"])
    return mathced_lines_idx
end


#loads_table

function _init_loads_df()
    return DataFrame(
        load_id=String[],
        source_id=String[],
        bus=String[],
        connections=Array{Int64,1}[],
        vm_nom=Float64[],
        pd_nom=Array{Float64,1}[],
        qd_nom=Array{Float64,1}[],
        configuration=String[],
        model=String[],
        dispatchable=String[],
        status=String[]
    )
end

function _push_loads_df!(loads_df, load_id, load)
    push!(loads_df, (
        load_id,
        get(load, "source_id", ""),
        get(load, "bus", ""),
        get(load, "connections", Int[]),
        get(load, "vm_nom", 0.0),
        get(load, "pd_nom", Float64[]),
        get(load, "qd_nom", Float64[]),
        string(get(load, "configuration", "")),
        string(get(load, "model", "")),
        string(get(load, "dispatchable", "")),
        string(get(load, "status", ""))
    ))
end

"""
    loads_table(eng::Dict{String, Any})

Generate a table summarizing the loads in the electrical network described by the dictionary `eng`.

# Arguments
- `eng::Dict{String, Any}`: A dictionary containing various components of the electrical network.

# Description
This function extracts the loads from the `eng` dictionary and creates a DataFrame with the load ID, status, bus, phases, kw, kvar, and kva. It then prints a formatted table of the loads.

# Example

```julia
using PowerModelsDistribution
using Pliers
eng= PowerModelsDistribution.parse_file("example.dss")
loads_table(eng)
```
"""
function loads_table(eng::Dict{String,Any})
    loads = haskey(eng, "load") ? eng["load"] : error("No loads found in the engineering data check the data model has \"load\" key")
    loads_df = _init_loads_df()

    for (load_id, load) in loads
        _push_loads_df!(loads_df, load_id, load)
    end
    header("Loads Table ($(nrow(loads_df)) loads)")
    pretty_table(sort!(loads_df))
    extra_keys(loads, ["load_id", "source_id", "status", "bus", "connections", "pd_nom", "qd_nom", "vm_nom", "model", "dispatchable", "configuration"])
end

""" 
    loads_table(eng::Dict{String,Any}, condition)
Generate and display a filtered table of loads from the given engineering data.

# Arguments
- `eng::Dict{String,Any}`: A dictionary containing engineering data, which must include a "load" key with load information.
- `condition`: A function that takes a load dictionary as input and returns a boolean indicating whether the load meets the filtering criteria.

# Description
This function extracts load information from the provided engineering data dictionary, applies the given condition to filter the loads, and then displays the filtered loads in a formatted table. Each load is augmented with its `load_id` before filtering. The table includes columns for `load_id`, `source_id`, `bus`, `connections`, `vm_nom`, `pd_nom`, `qd_nom`, `configuration`, `model`, `dispatchable`, and `status`.

# Example
- filter the loads by the value of `pd_nom`:
    ```julia
    loads_table(eng, load -> load["pd_nom"] > [0.33])
    ```
- filter the loads by the phase connectivity
    ```julia
    loads_table(eng, load -> load["connections"] == [1, 4])
    ```
"""
function loads_table(eng::Dict{String,Any}, condition::Function)
    loads = haskey(eng, "load") ? eng["load"] : error("No loads found in the engineering data check the data model has \"load\" key")
    # adding the load_id to the load dictionary so it can be filtered by name
    for (load_id, load) in loads
        load["load_id"] = load_id
    end

    loads_df = _init_loads_df()
    matched_loads_idx = []
    for (load_id, load) in loads
        if condition(load)
            _push_loads_df!(loads_df, load_id, load)
            push!(matched_loads_idx, load_id)
        end
    end
    header("Filtered Loads Table ($(nrow(loads_df)) loads)")
    pretty_table(sort!(loads_df))
    extra_keys(loads, ["load_id", "source_id", "status", "bus", "connections", "pd_nom", "qd_nom", "vm_nom", "model", "dispatchable", "configuration"])
    return matched_loads_idx
end



# linecode table


function linecodes_table(eng::Dict{String,Any})
    linecodes = haskey(eng, "linecode") ? eng["linecode"] : error("No linecodes found in the engineering data check the data model has \"linecode\" key")


    header("Linecodes Table ($(length(linecodes)) linecodes)")
    c = 0

    for (linecode_id, linecode) in linecodes
        c = c + 1
        sub_header("$(c) - Linecode: $(linecode_id)")

        for (key, matrix) in linecode
            sub_sub_header("$(linecode_id) - ($key) :")
            try
                _pretty_diag_matrix(matrix::Matrix)
            catch e
                error_text("parsing the $key it might not be a matrix: $(e)")
            end

        end
    end
    #pretty_table(sort!(linecodes_df))
    #extra_keys(linecodes, ["rs", "xs", "g_fr", "g_to", "b_fr", "b_to", "cm_ub"])
end

function linecode_table(eng::Dict{String,Any}, linecodes)
    haskey(eng, "linecode") ? eng["linecode"] : error("No linecodes found in the engineering data check the data model has \"linecode\" key")
    c = 0

    for linecode_id in linecodes
        c = c + 1
        sub_header("$(c) - Linecode: $(linecode_id)")

        for (key, matrix) in eng["linecode"][linecode_id]
            sub_sub_header("$(linecode_id) - ($key) :")
            try
                _pretty_diag_matrix(matrix::Matrix)
            catch e
                error_text("parsing the $key it might not be a matrix: $(e)")
            end

        end
    end
    #pretty_table(sort!(linecodes_df))
    #extra_keys(linecodes, ["rs", "xs", "g_fr", "g_to", "b_fr", "b_to", "cm_ub"])
end


## transformers_table

function _init_transformers_df()
    return DataFrame(
        transformer_id=String[],
        name=String[],
        source_id=String[],
        bus=Vector{String}[],
        connections=Vector{Vector{Int}}[],
        vm_nom=Vector{Float64}[],
        sm_nom=Vector{Float64}[],
        configuration=Vector{String}[],
        xsc=Vector{Float64}[],
        rw=Vector{Float64}[],
        noloadloss=Float64[],
        cmag=Float64[],
        tm_set=Vector{Vector{Float64}}[],
        tm_fix=Vector{Vector{Bool}}[],
        polarity=Vector{Int}[],
        status=String[]
    )
end

function _push_transformers_df!(transformers_df, transformer_id, transformer)
    push!(transformers_df, (
        transformer_id,
        get(transformer, "name", ""),
        get(transformer, "source_id", ""),
        get(transformer, "bus", String[]),
        get(transformer, "connections", Vector{Int}[]),
        get(transformer, "vm_nom", Float64[]),
        get(transformer, "sm_nom", Float64[]),
        [string(c) for c in get(transformer, "configuration", [])],
        get(transformer, "xsc", Float64[]),
        get(transformer, "rw", Float64[]),
        get(transformer, "noloadloss", 0.0),
        get(transformer, "cmag", 0.0),
        get(transformer, "tm_set", Vector{Float64}[]),
        get(transformer, "tm_fix", Vector{Bool}[]),
        get(transformer, "polarity", Int[]),
        string(get(transformer, "status", ""))
    ))
end


"""
    transformers_table(eng::Dict{String, Any})

Generate a table summarizing the transformers in the electrical network described by the dictionary `eng`.

# Arguments
- `eng::Dict{String, Any}`: A dictionary containing various components of the electrical network.

# Description
This function extracts the transformers from the `eng` dictionary and creates a DataFrame with the transformer ID, name, source ID, buses, connections, nominal voltage (vm_nom), nominal power (sm_nom), configuration, short-circuit reactance (xsc), winding resistance (rw), no-load loss, magnetizing current (cmag), tap settings (tm_set), fixed taps (tm_fix), polarity, and status. It then prints a formatted table of the transformers.

# Example
```julia
using PowerModelsDistribution
using Pliers 
eng= PowerModelsDistribution.parse_file("example.dss")
transformers_table(eng)
```
"""
function transformers_table(eng::Dict{String,Any})
    transformers = haskey(eng, "transformer") ? eng["transformer"] : error("No transformers found in the engineering data check the data model has \"transformer\" key")

    transformers_df = _init_transformers_df()

    for (transformer_id, transformer) in transformers
        _push_transformers_df!(transformers_df, transformer_id, transformer)
    end
    header("Transformers Table ($(nrow(transformers_df)) transformers)")
    pretty_table(sort!(transformers_df))

    extra_keys(transformers, ["transformer_id", "name", "source_id", "bus", "connections", "vm_nom", "sm_nom", "configuration", "xsc", "rw", "noloadloss", "cmag", "tm_set", "tm_fix", "tm_step", "polarity", "status"])
end

"""
    transformers_table(eng::Dict{String,Any}, condition)

Generate and display a filtered table of transformers from the given engineering data.

# Arguments
- `eng::Dict{String,Any}`: A dictionary containing engineering data, which must include a "transformer" key with transformer information.
- `condition`: A function that takes a transformer dictionary as input and returns a boolean indicating whether the transformer meets the filtering criteria.

# Description
This function extracts transformer information from the provided engineering data dictionary, applies the given condition to filter the transformers, and then displays the filtered transformers in a formatted table. Each transformer is augmented with its `transformer_id` before filtering.

# Example
```julia
using PowerModelsDistribution
using Pliers
eng= PowerModelsDistribution.parse_file("example.dss")
transformers_table(eng, transformer -> transformer["vm_nom"][1] > 10.0)
```
"""
function transformers_table(eng::Dict{String,Any}, condition::Function)
    transformers = haskey(eng, "transformer") ? eng["transformer"] : error("No transformers found in the engineering data check the data model has \"transformer\" key")
    # adding the transformer_id to the transformer dictionary so it can be filtered by name
    for (transformer_id, transformer) in transformers
        transformer["transformer_id"] = transformer_id
    end

    transformers_df = _init_transformers_df()
    matched_transformers_idx = []
    for (transformer_id, transformer) in transformers
        if condition(transformer)
            _push_transformers_df!(transformers_df, transformer_id, transformer)
            push!(matched_transformers_idx, transformer_id)
        end
    end
    header("Filtered Transformers Table ($(nrow(transformers_df)) transformers)")
    pretty_table(sort!(transformers_df))
    extra_keys(transformers, ["transformer_id", "name", "source_id", "bus", "connections", "vm_nom", "sm_nom", "configuration", "xsc", "rw", "noloadloss", "cmag", "tm_set", "tm_fix", "tm_step", "polarity", "status"])
    return matched_transformers_idx
end


## 

# Here I just want to add the ability for the above tables to be called even if the user passes a `.dss` file path instead of the parsed dictionary


macro define_table_from_file(name)
    quote
        """
        Processes a DSS file and calls the function with the parsed data.

        # Arguments
        - `dss::String`: The path to the DSS file to be processed.

        # Returns
        Whatever function returns when called with parsed data.

        """
        function $(esc(name))(dss::String)
            eng = PowerModelsDistribution.parse_file(dss)
            $(esc(name))(eng)
        end
    end
end

@define_table_from_file buses_table
@define_table_from_file lines_table
@define_table_from_file loads_table
@define_table_from_file linecodes_table
@define_table_from_file transformers_table



#=
░███     ░███    ░███    ░██████████░██     ░██    ░█████████  ░██████████ ░█████████    ░██████   ░█████████  ░██████████
░████   ░████   ░██░██       ░██    ░██     ░██    ░██     ░██ ░██         ░██     ░██  ░██   ░██  ░██     ░██     ░██    
░██░██ ░██░██  ░██  ░██      ░██    ░██     ░██    ░██     ░██ ░██         ░██     ░██ ░██     ░██ ░██     ░██     ░██    
░██ ░████ ░██ ░█████████     ░██    ░██████████    ░█████████  ░█████████  ░█████████  ░██     ░██ ░█████████      ░██    
░██  ░██  ░██ ░██    ░██     ░██    ░██     ░██    ░██   ░██   ░██         ░██         ░██     ░██ ░██   ░██       ░██    
░██       ░██ ░██    ░██     ░██    ░██     ░██    ░██    ░██  ░██         ░██          ░██   ░██  ░██    ░██      ░██    
░██       ░██ ░██    ░██     ░██    ░██     ░██    ░██     ░██ ░██████████ ░██           ░██████   ░██     ░██     ░██    



=#


# MATH report d=====(￣▽￣*)b

function math_report(math::Dict{String,Any}; detailed=false)
    print(UNDERLINE(BLUE_FG("Report for the ", BOLD("$(get(math, "data_model", nothing))"),
        " model of the network  ", BOLD("$(get(math, "name", nothing))"), " \n")))

    print(WHITE_FG("This network has:\n"))
    print("             $(BOLD("$(length(get(math, "gen", [])))")) generators, \n")
    print("             $(BOLD("$(length(get(math, "bus", [])))")) buses, \n")
    print("             $(BOLD("$(length(get(math, "branch", [])))")) branches, \n")
    print("             $(BOLD("$(length(get(math, "load", [])))")) loads, \n")
    print("             $(BOLD("$(length(get(math, "transformer", [])))")) transformers, \n")
    print("             $(BOLD("$(length(get(math, "time_series", [])))")) time series data points.\n")

    if detailed
        print("\n")
        math_buses_table(math)
        print("\n")
        math_branches_table(math)
        print("\n")
        math_loads_table(math)
        print("\n")
        math_gen_table(math)
        if haskey(math, "transformer") && length(math["transformer"]) > 0
            print("\n")
            math_transformers_table(math)
        end
    end

end
# MATH BUSES (. ❛ ᴗ ❛.)
function _initiate_math_buses_df()
    return DataFrame(
        index=Int[],
        bus_id=String[],
        bus_type=Int[],
        terminals=Vector{Int}[],
        grounded=Vector{Int}[],
        vm_pair_ub=Vector{Tuple{Any,Any,Real}}[],
        vm_pair_lb=Vector{Tuple{Any,Any,Real}}[],
        source_id=String[],
        vbase=Float64[],
        vmin=Vector{Float64}[],
        vmax=Vector{Float64}[],
    )
end

function _pushing_math_buses_df!(buses_df, bus_id, bus)
    push!(buses_df, (
        get(bus, "index", missing),
        bus_id,
        get(bus, "bus_type", missing),
        get(bus, "terminals", Int[]),
        get(bus, "grounded", Int[]),
        get(bus, "vm_pair_ub", Tuple{Any,Any,Real}[]),
        get(bus, "vm_pair_lb", Tuple{Any,Any,Real}[]),
        get(bus, "source_id", missing),
        get(bus, "vbase", missing),
        get(bus, "vmin", Float64[]),
        get(bus, "vmax", Float64[]),
    ))
end

function math_buses_table(math::Dict{String,Any})
    buses = haskey(math, "bus") ? math["bus"] : error("No buses found in the MATHEMATICAL data check the data model has \"bus\" key")
    buses_df = _initiate_math_buses_df()
    for (bus_id, bus) in buses
        _pushing_math_buses_df!(buses_df, bus_id, bus)
    end
    header("Buses Table ($(nrow(buses_df)) buses)")
    pretty_table(sort!(buses_df))
    extra_keys(buses, ["index", "bus_id", "bus_type", "terminals", "grounded", "vm_pair_ub", "vm_pair_lb", "source_id", "vbase", "vmin", "vmax"])
end

function math_buses_table(math::Dict{String,Any}, condition::Function)
    buses = haskey(math, "bus") ? math["bus"] : error("No buses found in the MATHEMATICAL data check the data model has \"bus\" key")
    buses_df = _initiate_math_buses_df()
    for (bus_id, bus) in buses
        if condition(bus)
            _pushing_math_buses_df!(buses_df, bus_id, bus)
        end
    end
    header("Buses Table ($(nrow(buses_df)) buses)")
    pretty_table(sort!(buses_df))
    extra_keys(buses, ["index", "bus_id", "bus_type", "terminals", "grounded", "vm_pair_ub", "vm_pair_lb", "source_id", "vbase", "vmin", "vmax"])
end

function math_bus_details(math, idx)
    labels = ["bus_id", "bus_type", "terminals", "grounded", "vm_pair_ub", "vm_pair_lb", "source_id", "vbase", "vmin", "vmax"]
    for id in idx
        bus = math["bus"][id]
        header("Bus $id")
        values = [get(bus, label, missing) for label in labels]
        table = DataFrame(labels=labels, values=values)
        pretty_table(table, header=["Name", "Value"])
    end
end

# MATH BRANCHES (. ❛ ᴗ ❛.) 

function _intitiate_math_branches_df()
    return DataFrame(
        index=Int[],
        source_id=String[],
        f_bus=Int[],
        t_bus=Int[],
        f_connections=Vector{Int}[],
        t_connections=Vector{Int}[],
        br_status=Int[],
        vbase=Float64[],
        br_r=Matrix{Float64}[],
        br_x=Matrix{Float64}[],
        g_fr=Matrix{Float64}[],
        b_fr=Matrix{Float64}[],
        g_to=Matrix{Float64}[],
        b_to=Matrix{Float64}[],
        rate_a=Vector{Real}[],
        rate_b=Vector{Real}[],
        rate_c=Vector{Real}[],
        c_rating_a=Vector{Real}[],
        c_rating_b=Vector{Real}[],
        c_rating_c=Vector{Real}[],
        angmin=Vector{Real}[],
        angmax=Vector{Real}[],
    )
end

function _pushing_math_branches_df!(branches_df, branch_id, branch)
    push!(branches_df, (
        get(branch, "index", missing),
        get(branch, "source_id", missing),
        get(branch, "f_bus", missing),
        get(branch, "t_bus", missing),
        get(branch, "f_connections", Int[]),
        get(branch, "t_connections", Int[]),
        get(branch, "br_status", missing),
        get(branch, "vbase", missing),
        get(branch, "br_r", Matrix{Float64}[]),
        get(branch, "br_x", Matrix{Float64}[]),
        get(branch, "g_fr", Matrix{Float64}[]),
        get(branch, "b_fr", Matrix{Float64}[]),
        get(branch, "g_to", Matrix{Float64}[]),
        get(branch, "b_to", Matrix{Float64}[]),
        get(branch, "rate_a", Real[]),
        get(branch, "rate_b", Real[]),
        get(branch, "rate_c", Real[]),
        get(branch, "c_rating_a", Real[]),
        get(branch, "c_rating_b", Real[]),
        get(branch, "c_rating_c", Real[]),
        get(branch, "angmin", Real[]),
        get(branch, "angmax", Real[]),
    ))

end


function math_branches_table(math::Dict{String,Any})
    branches = haskey(math, "branch") ? math["branch"] : error("No branches found in the MATHEMATICAL data check the data model has \"branch\" key")
    branches_df = _intitiate_math_branches_df()
    for (branch_id, branch) in branches
        _pushing_math_branches_df!(branches_df, branch_id, branch)
    end
    header("Branches Table ($(nrow(branches_df)) branches)")
    pretty_table(sort!(branches_df))
    extra_keys(branches, ["index", "name", "source_id", "f_bus", "t_bus", "f_connections", "t_connections", "br_status", "vbase", "br_r", "br_x", "g_fr", "b_fr", "g_to", "b_to", "rate_a", "rate_b", "rate_c", "c_rating_a", "c_rating_b", "c_rating_c", "angmin", "angmax"])
end

function math_branches_table(math::Dict{String,Any}, condition::Function)
    branches = haskey(math, "branch") ? math["branch"] : error("No branches found in the MATHEMATICAL data check the data model has \"branch\" key")
    branches_df = _intitiate_math_branches_df()
    matching_branch_idx = []
    for (branch_id, branch) in branches
        if condition(branch)
            _pushing_math_branches_df!(branches_df, branch_id, branch)
            push!(matching_branch_idx, branch_id)
        end
    end
    header("Branches Table ($(nrow(branches_df)) branches)")
    pretty_table(sort!(branches_df))
    extra_keys(branches, ["index", "name", "source_id", "f_bus", "t_bus", "f_connections", "t_connections", "br_status", "vbase", "br_r", "br_x", "g_fr", "b_fr", "g_to", "b_to", "rate_a", "rate_b", "rate_c", "c_rating_a", "c_rating_b", "c_rating_c", "angmin", "angmax"])
    return matching_branch_idx
end

function math_branch_details(math::Dict{String,Any}, branch_idx)
    branches = haskey(math, "branch") ? math["branch"] : error("No branches found in the MATHEMATICAL data check the data model has \"branch\" key")
    for id in branch_idx
        branch = branches[id]

        header("Branch $id")

        labels_1 = ["f_bus", "t_bus", "f_connections", "t_connections", "br_status", "vbase", "rate_a", "rate_b", "rate_c", "c_rating_a", "c_rating_b", "c_rating_c", "angmin", "angmax"]
        values_1 = [get(branch, label, missing) for label in labels_1]
        table_1 = DataFrame(labels=labels_1, values=values_1)
        pretty_table(table_1, header=["Name", "Value"])

        sub_header("br_r:  ")
        _pretty_diag_matrix(get(branch, "br_r", missing))
        sub_header("br_x:  ")
        _pretty_diag_matrix(get(branch, "br_x", missing))
        sub_header("g_fr:  ")
        _pretty_diag_matrix(get(branch, "g_fr", missing))
        sub_header("b_fr:  ")
        _pretty_diag_matrix(get(branch, "b_fr", missing))
        sub_header("g_to:  ")
        _pretty_diag_matrix(get(branch, "g_to", missing))
        sub_header("b_to:  ")
        _pretty_diag_matrix(get(branch, "b_to", missing))

    end
end



# MATH LOADS (^_~)  
function _initiate_math_loads_df()
    return DataFrame(
        index=Int[],
        source_id=String[],
        load_bus=Int[],
        connections=Vector{Int}[],
        configuration=String[],
        name=String[],
        pd=Vector{Real}[],
        qd=Vector{Real}[],
        status=Int[],
        vbase=Float64[],
        vnom_kv=Float64[],
        dispatchable=Int[],
    )
end

function _pushing_math_loads_df!(loads_df, load_id, load)
    push!(loads_df, (
        get(load, "index", missing),
        get(load, "source_id", missing),
        get(load, "load_bus", missing),
        get(load, "connections", Int[]),
        string(get(load, "configuration", missing)),
        get(load, "name", missing),
        get(load, "pd", Real[]),
        get(load, "qd", Real[]),
        get(load, "status", missing),
        get(load, "vbase", missing),
        get(load, "vnom_kv", missing),
        get(load, "dispatchable", missing),
    ))
end

function math_loads_table(math::Dict{String,Any})
    loads = haskey(math, "load") ? math["load"] : error("No loads found in the MATHEMATICAL data check the data model has \"load\" key")
    loads_df = _initiate_math_loads_df()
    for (load_id, load) in loads
        _pushing_math_loads_df!(loads_df, load_id, load)
    end
    header("Loads Table ($(nrow(loads_df)) loads)")
    pretty_table(sort!(loads_df))
    extra_keys(loads, ["index", "source_id", "load_bus", "connections", "configuration", "name", "status", "qd", "vbase", "vnom_kv", "dispatchable", "pd"])
end
function math_loads_table(math::Dict{String,Any}, condition::Function)
    loads = haskey(math, "load") ? math["load"] : error("No loads found in the MATHEMATICAL data check the data model has \"load\" key")
    loads_df = _initiate_math_loads_df()
    mathced_load_idx = []
    for (load_id, load) in loads
        if condition(load)
            _pushing_math_loads_df!(loads_df, load_id, load)
            push!(mathced_load_idx, load_id)
        end
    end
    header("Loads Table ($(nrow(loads_df)) loads)")
    pretty_table(sort!(loads_df))
    extra_keys(loads, ["index", "source_id", "load_bus", "connections", "configuration", "name", "status", "qd", "vbase", "vnom_kv", "dispatchable", "pd"])
    return mathced_load_idx
end



function math_load_details(math, idx)
    labels = ["load_bus", "connections", "configuration", "name", "pd", "qd", "status", "vbase", "vnom_kv", "dispatchable"]
    for id in idx
        load = math["load"][id]
        header("Load $id")
        values = [get(load, label, missing) for label in labels]
        table = DataFrame(labels=labels, values=values)
        pretty_table(table, header=["Name", "Value"])
    end
end


# MATH GEN (. ❛ ᴗ ❛.)


function _initiate_math_gen_df()
    return DataFrame(
        index=Int[],
        gen_bus=Int[],
        source_id=String[],
        connections=Vector{Int}[],
        configuration=String[],
        name=String[],
        pg=Vector{Real}[],
        qg=Vector{Real}[],
        pmax=Vector{Real}[],
        pmin=Vector{Real}[],
        qmax=Vector{Real}[],
        qmin=Vector{Real}[],
        vg=Vector{Real}[],
        vbase=Float64[],
        gen_status=Int[],
        model=Int[],
        cost=Vector{Real}[],
        ncost=Int[],
        control_mode=Int[],
        startup=Float64[],
        shutdown=Float64[],
    )
end

function _pushing_math_gen_df!(gen_df, gen_id, gen)
    push!(gen_df, (
        get(gen, "index", missing),
        get(gen, "gen_bus", missing),
        get(gen, "source_id", missing),
        get(gen, "connections", Int[]),
        string(get(gen, "configuration", missing)),
        get(gen, "name", missing),
        get(gen, "pg", Real[]),
        get(gen, "qg", Real[]),
        get(gen, "pmax", Real[]),
        get(gen, "pmin", Real[]),
        get(gen, "qmax", Real[]),
        get(gen, "qmin", Real[]),
        get(gen, "vg", Real[]),
        get(gen, "vbase", missing),
        get(gen, "gen_status", missing),
        get(gen, "model", missing),
        get(gen, "cost", Vector{Real}[]),
        get(gen, "ncost", missing),
        get(gen, "control_mode", missing),
        get(gen, "startup", Real[]),
        get(gen, "shutdown", Real[]),
    ))
end

function math_gen_table(math::Dict{String,Any})
    gens = haskey(math, "gen") ? math["gen"] : error("No generators found in the MATHEMATICAL data check the data model has \"gen\" key")
    gens_df = _initiate_math_gen_df()
    for (gen_id, gen) in gens
        _pushing_math_gen_df!(gens_df, gen_id, gen)
    end
    header("Generators Table ($(nrow(gens_df)) generators)")
    pretty_table(sort!(gens_df))
    extra_keys(gens, ["index", "source_id", "gen_bus", "connections", "configuration", "name", "pg", "qg", "pmax", "pmin", "qmax", "qmin", "vg", "vbase", "gen_status", "model", "cost", "ncost", "control_mode", "startup", "shutdown"])
end

function math_gen_table(math::Dict{String,Any}, condition::Function)
    gens = haskey(math, "gen") ? math["gen"] : error("No generators found in the MATHEMATICAL data check the data model has \"gen\" key")
    gens_df = _initiate_math_gen_df()
    mathced_gen_idx = []
    for (gen_id, gen) in gens
        if condition(gen)
            _pushing_math_gen_df!(gens_df, gen_id, gen)
            push!(mathced_gen_idx, gen_id)
        end
    end
    header("Generators Table ($(nrow(gens_df)) generators)")
    pretty_table(sort!(gens_df))
    extra_keys(gens, ["index", "source_id", "gen_bus", "connections", "configuration", "name", "pg", "qg", "pmax", "pmin", "qmax", "qmin", "vg", "vbase", "gen_status", "model", "cost", "ncost", "control_mode", "startup", "shutdown"])
    return mathced_gen_idx
end

function math_gen_details(math, idx)

    labels = ["gen_bus", "source_id", "connections", "configuration", "name", "pg", "qg", "pmax", "pmin", "qmax", "qmin", "vg", "vbase", "gen_status", "model", "cost", "ncost", "control_mode", "startup", "shutdown"]
    for id in idx
        gen = math["gen"][id]
        header("Generator $id")
        values = [get(gen, label, missing) for label in labels]
        table = DataFrame(labels=labels, values=values)
        pretty_table(table, header=["Name", "Value"])
    end
end


# MATH TRANSFORMERS (⚡)

function _initiate_math_transformers_df()
    return DataFrame(
        index=Int[],
        name=String[],
        source_id=String[],
        f_bus=Int[],
        t_bus=Int[],
        f_connections=Vector{Int}[],
        t_connections=Vector{Int}[],
        configuration=String[],
        tm_nom=Float64[],
        tm_set=Vector{Float64}[],
        tm_fix=Vector{Bool}[],
        tm_step=Vector{Float64}[],
        f_vbase=Float64[],
        t_vbase=Float64[],
        polarity=Int[],
        sm_ub=Float64[],
        cm_ub=Float64[],
        status=Int[],
    )
end

function _pushing_math_transformers_df!(transformers_df, transformer_id, transformer)
    push!(transformers_df, (
        get(transformer, "index", missing),
        get(transformer, "name", ""),
        get(transformer, "source_id", ""),
        get(transformer, "f_bus", missing),
        get(transformer, "t_bus", missing),
        get(transformer, "f_connections", Int[]),
        get(transformer, "t_connections", Int[]),
        string(get(transformer, "configuration", "")),
        get(transformer, "tm_nom", missing),
        get(transformer, "tm_set", Float64[]),
        get(transformer, "tm_fix", Bool[]),
        get(transformer, "tm_step", Float64[]),
        get(transformer, "f_vbase", missing),
        get(transformer, "t_vbase", missing),
        get(transformer, "polarity", missing),
        get(transformer, "sm_ub", missing),
        get(transformer, "cm_ub", missing),
        get(transformer, "status", missing),
    ))
end

"""
    math_transformers_table(math::Dict{String, Any})

Generate a table summarizing the transformers in the mathematical model described by the dictionary `math`.

# Arguments
- `math::Dict{String, Any}`: A dictionary containing various components of the mathematical model.

# Description
This function extracts the transformers from the `math` dictionary and creates a DataFrame with the transformer index, name, source ID, from/to buses, connections, configuration, tap settings, voltage bases, polarity, power/current limits, and status. It then prints a formatted table of the transformers.

# Example
```julia
using PowerModelsDistribution
using Pliers 
eng = PowerModelsDistribution.parse_file("example.dss")
math = PowerModelsDistribution.transform_data_model(eng)
math_transformers_table(math)
```
"""
function math_transformers_table(math::Dict{String,Any})
    transformers = haskey(math, "transformer") ? math["transformer"] : error("No transformers found in the MATHEMATICAL data check the data model has \"transformer\" key")
    transformers_df = _initiate_math_transformers_df()
    for (transformer_id, transformer) in transformers
        _pushing_math_transformers_df!(transformers_df, transformer_id, transformer)
    end
    header("Transformers Table ($(nrow(transformers_df)) transformers)")
    pretty_table(sort!(transformers_df))
    extra_keys(transformers, ["index", "name", "source_id", "f_bus", "t_bus", "f_connections", "t_connections", "configuration", "tm_nom", "tm_set", "tm_fix", "tm_step", "f_vbase", "t_vbase", "polarity", "sm_ub", "cm_ub", "status"])
end

"""
    math_transformers_table(math::Dict{String, Any}, condition::Function)

Generate and display a filtered table of transformers from the given mathematical model data.

# Arguments
- `math::Dict{String,Any}`: A dictionary containing mathematical model data, which must include a "transformer" key with transformer information.
- `condition`: A function that takes a transformer dictionary as input and returns a boolean indicating whether the transformer meets the filtering criteria.

# Description
This function extracts transformer information from the provided mathematical model dictionary, applies the given condition to filter the transformers, and then displays the filtered transformers in a formatted table.

# Example
```julia
using PowerModelsDistribution
using Pliers
eng = PowerModelsDistribution.parse_file("example.dss")
math = PowerModelsDistribution.transform_data_model(eng)
math_transformers_table(math, transformer -> transformer["f_bus"] == 1)
```
"""
function math_transformers_table(math::Dict{String,Any}, condition::Function)
    transformers = haskey(math, "transformer") ? math["transformer"] : error("No transformers found in the MATHEMATICAL data check the data model has \"transformer\" key")
    transformers_df = _initiate_math_transformers_df()
    matched_transformer_idx = []
    for (transformer_id, transformer) in transformers
        if condition(transformer)
            _pushing_math_transformers_df!(transformers_df, transformer_id, transformer)
            push!(matched_transformer_idx, transformer_id)
        end
    end
    header("Transformers Table ($(nrow(transformers_df)) transformers)")
    pretty_table(sort!(transformers_df))
    extra_keys(transformers, ["index", "name", "source_id", "f_bus", "t_bus", "f_connections", "t_connections", "configuration", "tm_nom", "tm_set", "tm_fix", "tm_step", "f_vbase", "t_vbase", "polarity", "sm_ub", "cm_ub", "status"])
    return matched_transformer_idx
end

function math_transformer_details(math, idx)
    labels = ["name", "source_id", "f_bus", "t_bus", "f_connections", "t_connections", "configuration", "tm_nom", "tm_set", "tm_fix", "tm_step", "f_vbase", "t_vbase", "polarity", "sm_ub", "cm_ub", "status"]
    for id in idx
        transformer = math["transformer"][id]
        header("Transformer $id")
        values = [get(transformer, label, missing) for label in labels]
        table = DataFrame(labels=labels, values=values)
        pretty_table(table, header=["Name", "Value"])
    end
end










#=
░██████████ ░██    ░██ ░█████████    ░██████   ░█████████  ░██████████
░██          ░██  ░██  ░██     ░██  ░██   ░██  ░██     ░██     ░██    
░██           ░██░██   ░██     ░██ ░██     ░██ ░██     ░██     ░██    
░█████████     ░███    ░█████████  ░██     ░██ ░█████████      ░██    
░██           ░██░██   ░██         ░██     ░██ ░██   ░██       ░██    
░██          ░██  ░██  ░██          ░██   ░██  ░██    ░██      ░██    
░██████████ ░██    ░██ ░██           ░██████   ░██     ░██     ░██    



=#

# Re-export PMD utility functions from parent module
export fluff_bus_voltages!
export solution_dictify_buses!, solution_dictify_branches!, solution_dictify_loads!, dictify_solution!
export calc_bases_from_dict, add_vmn_p_q
export kron_reduce_impedance, get_sequence_components
export show_example, show_transformer_math_components

# Re-export result explorer functions
export pf_results, bus_res, branch_viz

# Re-export eng explorer functions
export eng_report, buses_table, lines_table, loads_table, linecodes_table, linecode_table, transformers_table

# Re-export math explorer functions
export math_report, math_buses_table, math_branches_table, math_branch_details
export math_loads_table, math_gen_table, math_gen_details, math_load_details
export math_bus_details, math_transformers_table, math_transformer_details

# Re-export network graph functions
export get_graph_node, get_graph_edge, create_network_graph

end # module PMDUtils
