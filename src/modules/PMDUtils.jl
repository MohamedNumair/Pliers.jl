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



# pretty terminal packages
using Crayons
using Crayons.Box
using PrettyTables

# plotting packages
using Makie
# using MakieCore
using CairoMakie
using WGLMakie
if Sys.iswindows()
    # using GLMakie
end

# using Tyler
# using Tyler.TileProviders
# using Tyler.MapTiles
# using Tyler.Extents

# data analysis packages
using Statistics
using LinearAlgebra

using DataFrames
using CSV

include("../core/styles.jl")
include("../core/utils.jl")


_N_IDX = 4 # representing the neutral index in terminals 

# using Dates
# using StatsPlots
# using StatsBase
# using Distributions
# using Random  

# Power Distribution Tools
#using PowerModelsDistribution
#using PowerModelsDistributionStateEstimation
#using Ipopt


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

    
    if detailed == true
        _eng_extra_details_table(eng)
        buses_table(eng)
        lines_table(eng)
        loads_table(eng)
        if haskey(eng, "transformer") && length(eng["transformer"]) > 0
            transformers_table(eng)
        end
    end

end

function _eng_extra_details_table(eng::Dict{String,Any})
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
        #"Number of nodes" n_nodes;
        #"Number of branches" n_lines;
        #"Number of loads" n_loads;
        #"Number Transformers" n_tfs;
        "Total resistance (Ω)" round(tot_r, digits=4);
        "Total reactance (Ω)" round(tot_x, digits=4);
        "Total MVA of transformers" round(tot_mva_val, digits=4);
        "Total length of lines (m)" round(tot_len, digits=4);
        #"Area (sq km)"              round(area, digits=4);
        "Number of line types used" n_linecodes;
        "Total impedance (Ω)" round(z_tot, digits=4);
        "R/X ratio" round(rx, digits=4);
        "Length/load" round(len_load, digits=4);
        "Z/ load" round(z_cons, digits=4)
        #"GIS information"           gis_info;
        #"Max current capacity"      round(max_amps, digits=4);
    ]
    #display(data)
    #println("\"Feeder name\", \"Number of nodes\", \"Number of branches\", \"Number of loads\", \"Number Transformers\", \"Total resistance\", \"Total reactance\", \"Total MVA of transformers\", \"Total length of lines (m)\", \"Number of line types used\", \"Total impedance\", \"R/X ratio\", \"Length/load\", \"Z/ consumer\"")
    #println("$(eng["name"]),$(n_nodes), $(n_lines), $(n_loads), $(n_tfs), $(round(tot_r, digits=4)), $(round(tot_x, digits=4)), $(round(tot_mva_val, digits=4)), $(round(tot_len, digits=4)), $(n_linecodes), $(round(z_tot, digits=4)), $(round(rx, digits=4)), $(round(len_load, digits=4)), $(round(z_cons, digits=4))")
    #return data
    header("Meta Information Table")
    pretty_table(DataFrame(Attributes=data[:, 1], Values=data[:, 2]); alignment=[:l, :r])
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
░█████████                               ░██████████                       ░██                                         
░██     ░██                              ░██                               ░██                                         
░██     ░██  ░███████   ░███████         ░██         ░██    ░██ ░████████  ░██  ░███████  ░██░████  ░███████  ░██░████ 
░█████████  ░██    ░██ ░██               ░█████████   ░██  ░██  ░██    ░██ ░██ ░██    ░██ ░███     ░██    ░██ ░███     
░██   ░██   ░█████████  ░███████         ░██           ░█████   ░██    ░██ ░██ ░██    ░██ ░██      ░█████████ ░██      
░██    ░██  ░██               ░██        ░██          ░██  ░██  ░███   ░██ ░██ ░██    ░██ ░██      ░██        ░██      
░██     ░██  ░███████   ░███████  ░██    ░██████████ ░██    ░██ ░██░█████  ░██  ░███████  ░██       ░███████  ░██      
                                                                ░██                                                    
                                                                ░██                                                    
                                                                                                                       
=#
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
function calc_bases_from_dict(data::Dict{String,Any}; return_dict = true)   
    

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

    @debug "Creating bus phasor plot for bus ID: $bus_id"
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


function bus_phasor(bus_data::Dict{String, Any};
                    makie_backend=WGLMakie, fig_size=(800, 800), colors=[:darkred, :darkgreen, :darkblue, :black], kwargs...)
                    
    makie_backend.activate!()
    haskey(bus_data, "voltage") || error("The bus_data dictionary must contain a 'bus' dictionary with 'voltage' entry of the complex values. see also `dictify_solution!` function")
    f = Figure(size=fig_size)
    degree_ticks = 0:30:330
    radian_ticks = deg2rad.(degree_ticks)
    tick_labels = ["$(d)°" for d in degree_ticks]
    ax = PolarAxis(f[1,1],
                   title = "Bus Voltage Phasors",
                   thetaticks = (radian_ticks, tick_labels),
                   kwargs...
                  )
    for (i, color) in enumerate(colors)
        phase_key = string(i)
        if haskey(bus_data["voltage"], phase_key)
            V_complex = bus_data["voltage"][phase_key]
            Vm = abs(V_complex)
            θ = angle(V_complex)
            lines!(ax, [0, θ], [0, Vm], color=color, linewidth=2, linestyle=:solid; kwargs...)
        end
    end
    display(f)
    return f, ax
end


function bus_phasor!(ax::PolarAxis,bus_data::Dict{String, Any};
                     colors=[:darkred, :darkgreen, :darkblue, :black], kwargs...)
    haskey(bus_data, "voltage") || error("The bus_data dictionary must contain a 'bus' dictionary with 'voltage' entry of the complex values. see also `dictify_solution!` function")
    
    for (i, color) in enumerate(colors)
        phase_key = string(i)
        if haskey(bus_data["voltage"], phase_key)
            V_complex = bus_data["voltage"][phase_key]
            Vm = abs(V_complex)
            θ = angle(V_complex)
            lines!(ax, [0, θ], [0, Vm], color=color ; kwargs...)
        end
    end
    
    return ax
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
@recipe(VPhasor, value) do scene
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

#=

░███     ░███ ░██                             ░██     ░██    ░██    ░██░██            
░████   ░████                                 ░██     ░██    ░██       ░██            
░██░██ ░██░██ ░██ ░███████   ░███████         ░██     ░██ ░████████ ░██░██  ░███████  
░██ ░████ ░██ ░██░██        ░██    ░██        ░██     ░██    ░██    ░██░██ ░██        
░██  ░██  ░██ ░██ ░███████  ░██               ░██     ░██    ░██    ░██░██  ░███████  
░██       ░██ ░██       ░██ ░██    ░██         ░██   ░██     ░██    ░██░██        ░██ 
░██       ░██ ░██ ░███████   ░███████  ░██      ░██████       ░████ ░██░██  ░███████  
                                                                                      
=#                                                                                    
                                                                                      


"""
    remove_virtual_bus!(math_b)

Remove virtual bus, branch, and generator from the given `math` dictionary.

# Arguments
- `math_b::Dict`: A dictionary containing information about buses, branches, and generators.

# Example
"""
function remove_virtual_bus!(math)
    virtual_bus = findfirst(bus -> contains(bus["name"], "virtual_bus.voltage_source.source"), math["bus"])
    virtual_branch = findfirst(branch -> contains(branch["name"], "_virtual_branch.voltage_source.source"), math["branch"])
    virtual_gen = findfirst(gen -> contains(gen["name"], "_virtual_gen"), math["gen"])
    
    r_new = math["branch"][virtual_branch]["t_bus"]
    
    # math["bus"][string(r_new)]["va"] = math["bus"][string(virtual_bus)]["va"]
    # math["bus"][string(r_new)]["vm"] = math["bus"][string(virtual_bus)]["vm"]
    # math["bus"][string(r_new)]["vmin"] = math["bus"][string(virtual_bus)]["vmin"]
    # math["bus"][string(r_new)]["vmax"] = math["bus"][string(virtual_bus)]["vmax"]
    math["bus"][string(r_new)]["bus_type"] = 3
    math["gen"][string(virtual_gen)]["gen_bus"] = r_new

    delete!(math["bus"], string(virtual_bus))
    delete!(math["branch"], string(virtual_branch))
    

    #delete!(math_b["gen"], string(virtual_gen))

    return virtual_bus, r_new
end


function _is_eng(data::Dict{String, Any})
    if haskey(data, "data_model")
        if string(data["data_model"]) == "ENGINEERING"
            return true
        elseif string(data["data_model"]) == "MATHEMATICAL"
            return false
        else
            error("Invalid data model found in the provided model dictionary. it has to be either ENGINEERING or MATHEMATICAL")
        end
    else
        error("No data model found in the provided model dictionary.")
    end
end




### Operation on Results Dictionary ###
"""
    fluff_bus_voltages(PF_Res::Dict{String, Any})

This function is used to get all possible voltage forms from the output of the OPF voltage variables, so it processes the bus voltages in the given Power Flow Results dictionary `PF_Res`. For each bus in `PF_Res["bus"]`, it checks if the bus has real (`vr`) and imaginary (`vi`) voltage components. If both components are present, it calculates the complex voltage `V`, the voltage magnitude `vm`, and the voltage angle `va` for the bus. A warning is issued indicating that only `vr` and `vi` can be processed.

# Arguments
- `PF_Res::Dict{String, Any}`: A dictionary containing the Power Flow Results solution, including bus voltage information.
!!! You need to pass the ["solution"] key of the Power Flow Results dictionary.

# Modifies
- Adds the following keys to each bus dictionary if `vr` and `vi` are present:
  - `V`: The complex voltage calculated as `vr + vi*im`.
  - `vm`: The magnitude of the complex voltage.
  - `va`: The angle of the complex voltage.
"""
function fluff_bus_voltages!(data::Dict{String, Any}) 
    data = haskey(data,"solution") ? data["solution"] : data

    for (_, bus) in data["bus"]
        if  haskey(bus, "vr") && haskey(bus, "vi") 
            bus["V"] = bus["vr"] .+ bus["vi"]*im
            bus["vm"] = abs.(bus["V"])
            bus["va"] = angle.(bus["V"])
        elseif haskey(bus, "va") && haskey(bus, "vm") #TODO: check if when you apply the `solution_make_si` does the angles become in degrees ?
            bus["V"] = bus["vm"] .* exp.(im.*bus["va"])
            bus["vr"] = real.(bus["V"])
            bus["vi"] = imag.(bus["V"])
        else 
            warning_text("Can't fluff $(keys(bus)). Now I can just fluff the (`vr`, `vi`) and (`va`, `vm`) pairs.")
        end 
    end

end 


function solution_dictify_buses!(pf_sol::Dict{String, Any}, math::Dict{String, Any};  formulation = "IVR")
    fluff_bus_voltages!(pf_sol)
    pf_sol = haskey(pf_sol,"solution") ? pf_sol["solution"] : pf_sol

    for (b, bus) in math["bus"]
        terminals = math["bus"][b]["terminals"]
        # write a dictionary where the key is the terminal number and the value is the voltage at that terminal
        pf_sol["bus"][b]["voltage"] = Dict(string(term) =>  pf_sol["bus"][b]["V"][i] for (i, term) in enumerate(terminals))
    end 

end

function solution_dictify_loads!(pf_sol::Dict{String, Any}, math::Dict{String, Any}; formulation = "IVR") 
    # The idea is to create the complex current and power for each load and store it in the load dictionary under the key "current" and "power" respectively.
    # The current is calculated as `I = P + Q*im` and the power is calculated as `S = P + Q*im`  
    for (l, load) in pf_sol["load"]
        terminals = math["load"][l]["connections"]

        if !(string(math["load"][l]["configuration"]) == "DELTA")
            
            # Create current dictionary based on available keys
            if haskey(load, "crd_bus") && haskey(load, "cid_bus")
                load["current_bus"] = Dict(string(term) => load["crd_bus"][i] + load["cid_bus"][i]*im for (i, term) in enumerate(terminals))
                load["power_bus"] = Dict(string(term) => pf_sol["bus"][string(math["load"][l]["load_bus"])]["voltage"][string(term)] * load["current_bus"][string(term)] for (i, term) in enumerate(terminals))
            elseif haskey(load, "crd") && haskey(load, "cid")
                load["current"] = Dict(string(term) => load["crd"][i] + load["cid"][i]*im for (i, term) in enumerate(terminals))
                load["power"] = Dict(string(term) => pf_sol["bus"][string(math["load"][l]["load_bus"])]["voltage"][string(term)] * load["current"][string(term)] for (i, term) in enumerate(terminals))
            end

            # Create power dictionary based on available keys
            if haskey(load, "pd_bus") && haskey(load, "qd_bus")
                load["power"] = Dict(string(term) => load["pd_bus"][i] + load["qd_bus"][i]*im for (i, term) in enumerate(terminals))
                @debug "added bus powers"
            elseif haskey(load, "pd") && haskey(load, "qd")
                load["power"] = Dict(string(term) => load["pd"][i] + load["qd"][i]*im for (i, term) in enumerate(terminals))
                @debug "added load powers"
            end  
        else
            terminals = setdiff(terminals, _N_IDX)
                # Create current dictionary based on available keys
                if haskey(load, "crd_bus") && haskey(load, "cid_bus")
                    load["current_bus"] = Dict(string(term) => load["crd_bus"][i] + load["cid_bus"][i]*im for (i, term) in enumerate(terminals))
                end

                if haskey(load, "crd") && haskey(load, "cid")
                    load["current"] = Dict(string(term) => load["crd"][i] + load["cid"][i]*im for (i, term) in enumerate(terminals))
                end
    
                # Create power dictionary based on available keys
                # if haskey(load, "pd_bus") && haskey(load, "qd_bus")
                #     load["power_bus"] = Dict(string(term) => load["pd_bus"][i] + load["qd_bus"][i]*im for (i, term) in enumerate(terminals))
                #     @debug "DELTA: added bus powers"
                # end
                if haskey(load, "pd") && haskey(load, "qd")
                    load["power"] = Dict(string(term) => load["pd"][i] + load["qd"][i]*im for (i, term) in enumerate(terminals))
                end  
        end

    end
end 

function solution_dictify_branches!(pf_sol::Dict{String, Any}, math::Dict{String, Any}; formulation = "IVR")
    # The idea is to create the complex current and power for each branch and store it in the branch dictionary under the key "current" and "power" respectively.
    # The current is calculated as `I = P + Q*im` and the power is calculated as `S = P + Q*im`  
    for (b, branch) in math["branch"]
        f_terminals = math["branch"][b]["f_connections"]
        t_terminals = math["branch"][b]["t_connections"]
        # write a dictionary where the key is the terminal number and the value is the current at that terminal
        branch = pf_sol["branch"][b]
        branch["power_from"] = Dict(string(term) => branch["pf"][i] + branch["qf"][i]*im for (i, term) in enumerate(f_terminals))
        branch["power_to"] = Dict(string(term) => branch["pt"][i] + branch["qt"][i]*im for (i, term) in enumerate(t_terminals))
        #TODO: I can make it check the branch f_bus and t_bus and get their voltages and then calcualte the currents from the powers.

        haskey(branch, "csr_to") ? branch["shunt_current_to"] = Dict(string(term) => branch["csr_to"][i] + branch["csi_to"][i]*im for (i, term) in enumerate(t_terminals)) : nothing
        haskey(branch, "cr_to") ? branch["current_to"] = Dict(string(term) => branch["cr_to"][i] + branch["ci_to"][i]*im for (i, term) in enumerate(t_terminals)) : nothing
        haskey(branch, "csr_fr") ? branch["shunt_current_from"] = Dict(string(term) => branch["csr_fr"][i] + branch["csi_fr"][i]*im for (i, term) in enumerate(f_terminals)) : nothing
        haskey(branch, "cr_fr") ? branch["current_from"] = Dict(string(term) => branch["cr_fr"][i] + branch["ci_fr"][i]*im for (i, term) in enumerate(f_terminals))    : nothing
    
    end
end

function solution_dictify_gens!(pf_sol::Dict{String, Any}, math::Dict{String, Any}; formulation = "IVR")
    # The idea is to create the complex current and power for each generator and store it in the generator dictionary under the key "current" and "power" respectively.
    # The current is calculated as `I = P + Q*im` and the power is calculated as `S = P + Q*im`  
    for (g, gen) in pf_sol["gen"]
        terminals = setdiff(math["gen"][g]["connections"], [4]) 
        # write a dictionary where the key is the terminal number and the value is the current at that terminal
        haskey(gen, "crg") || haskey(gen, "cig") ?  gen["current"] = Dict(string(term) => gen["crg"][i] + gen["cig"][i]*im for (i, term) in enumerate(terminals)) : nothing
        haskey(gen, "pg") || haskey(gen, "pg_bus") ? gen["power"] = Dict(string(term) => gen["pg"][i] + gen["qg"][i]*im for (i, term) in enumerate(terminals)) : nothing
    end
end

"""
    dictify_solution!(pf_sol::Dict{String, Any}, math::Dict{String, Any}; formulation = "IVR")

Transforms and organizes the solution data from a power flow computation into a structured dictionary format.

# Arguments
- `pf_sol::Dict{String, Any}`: A dictionary containing the power flow solution data.
- `math::Dict{String, Any}`: A dictionary containing the mathematical model data.
- `formulation::String` (optional): Specifies the formulation type to be used. Defaults to `"IVR"`.

# Description
This function modifies the `pf_sol` dictionary in-place by calling helper functions to process and structure
data for buses, loads, branches, and generators. Each helper function is responsible for handling a specific
component of the power flow solution.

# Notes
- The function assumes that the helper functions `solution_dictify_buses!`, `solution_dictify_loads!`,
  `solution_dictify_branches!`, and `solution_dictify_gens!` are defined and properly handle their respective
  components.
- The `formulation` parameter allows customization of the solution processing based on the formulation type.
"""
function dictify_solution!(pf_sol::Dict{String, Any}, math::Dict{String, Any}; formulation = "IVR")
    pf_sol = haskey(pf_sol,"solution") ? pf_sol["solution"] : pf_sol
    solution_dictify_buses!(pf_sol, math; formulation = formulation)
    solution_dictify_loads!(pf_sol, math; formulation = formulation)
    solution_dictify_branches!(pf_sol, math; formulation = formulation)
    solution_dictify_gens!(pf_sol, math; formulation = formulation)
end

"""
    _separate_phase_neutral_voltages(pf_sol, bus_index)

Separate the phase and neutral voltages from the power flow solution for a given bus.

# Arguments
- `pf_sol::Dict`: The power flow solution dictionary containing voltage information.
- `bus_index::Int`: The index of the bus for which to separate the voltages.

# Returns
- `phase_voltage::Vector{ComplexF64}`: A vector containing the phase voltages (up to 3 phases).
- `neutral_voltage::ComplexF64`: The neutral voltage. If the neutral voltage is not present, returns 0 + 0im "assuming grounded".

"""
function _separate_phase_neutral_voltages(pf_sol, bus_index)
    phase_voltage = ComplexF64[]
    for i in 1:3 
        if haskey(pf_sol["bus"][bus_index]["voltage"], string(i))
            push!(phase_voltage, pf_sol["bus"][bus_index]["voltage"][string(i)])
        end
    end

    neutral_voltage = ComplexF64[]

    if haskey(pf_sol["bus"][bus_index]["voltage"], string(_N_IDX))  
        neutral_voltage = pf_sol["bus"][bus_index]["voltage"][string(_N_IDX)]
    else 
        neutral_voltage = 0 + 0im # assuming grounded
    end

    return phase_voltage, neutral_voltage
end



function _get_vmn(pf_sol, math, math_meas)
    for (b, bus) in pf_sol["bus"]
        pvs, vn = _separate_phase_neutral_voltages(pf_sol, b)
        vmn = abs.(pvs .- vn)
        bus["vmn"] = vmn
        math_meas["bus"][b]["terminals"] = setdiff(math["bus"][b]["terminals"], _N_IDX)
    end
end

function _get_pd_qd(pf_sol, math, math_meas)
    for (l, _) in pf_sol["load"]
        math_meas["load"][l]["connections"] = setdiff(math["load"][l]["connections"], _N_IDX)
    end
end

"""
    add_vmn_p_q(math, pf_sol) -> math_meas

Adds voltage magnitude (VMN), active power (P), and reactive power (Q) measurements 
to the given mathematical model.

# Arguments
- `math`: The mathematical model to which the measurements will be added.
- `pf_sol`: The power flow solution containing the necessary data for VMN, P, and Q.

# Returns
- `math_meas`: A deep copy of the input `math` with the added VMN, P, and Q measurements.

# Notes
This function internally calls `_get_vmn` to add voltage magnitude measurements 
and `_get_pd_qd` to add active and reactive power measurements.
"""
function add_vmn_p_q(math, pf_sol)
    math_meas = deepcopy(math)
    _get_vmn(pf_sol, math, math_meas)
    _get_pd_qd(pf_sol, math, math_meas)
    _add_delta_readings(pf_sol, math, math_meas)
    return math_meas
end


function _add_delta_readings(pf_sol, math, math_meas)

    for (l, load) in math["load"]
        if load["configuration"] == DELTA
            pf_sol["load"][l]["ptot"] = [sum(real(value) for (key, value) in pf_sol["load"][l]["power"])]
            pf_sol["load"][l]["qtot"] = [sum(imag(value) for (key, value) in pf_sol["load"][l]["power"])]
            
            #pf_sol["load"][l]["ptot"] = [sum(real(value) for (key, value) in pf_sol["load"][l]["power_bus"])]
            #pf_sol["load"][l]["qtot"] = [sum(imag(value) for (key, value) in pf_sol["load"][l]["power_bus"])]


            load_bus = pf_sol["bus"][string(load["load_bus"])]
            load_bus["vll"] = sqrt.([ (load_bus["vr"][x] - load_bus["vr"][y] )^2 + (load_bus["vi"][x] - load_bus["vi"][y] )^2 for (x,y) in [(1,2), (2,3), (3,1)]])
        end
    end
return math_meas
end

function write_delta_readings(pf_sol, math, math_meas)

    # for (b,bus) in pf_sol["bus"]
    #     bus["vmn2"] = []
    #     for t in setdiff(math["bus"][b]["terminals"],[4])
    #         push!(bus["vmn2"], (bus["vr"][t] - bus["vr"][4])^2 + (bus["vi"][t] - bus["vi"][4])^2)
    #     end
    #     bus["vmn"] = sqrt.(bus["vmn2"])
    # math_meas["bus"][b]["terminals"] = setdiff(math["bus"][b]["terminals"], 4)
    # end 

    for (l, load) in math["load"]
        if load["configuration"] == DELTA
            pf_sol["load"][l]["ptot"] = [sum(pf_sol["load"][l]["pd"])]
            pf_sol["load"][l]["qtot"] = [sum(pf_sol["load"][l]["qd"])]

            load_bus = pf_sol["bus"][string(load["load_bus"])]
            load_bus["vll"] = sqrt.([ (load_bus["vr"][x] - load_bus["vr"][y] )^2 + (load_bus["vi"][x] - load_bus["vi"][y] )^2 for (x,y) in [(1,2), (2,3), (3,1)]])
        end
    end
return math_meas
end


function ptot_qtot_from_loads(math, pf_sol)
    ptot = 0.0
    qtot = 0.0
    for (l, load) in pf_sol["load"]
        ptot += sum(real.(values(load["power"])))
        qtot += sum(imag.(values(load["power"])))
    end
    return ptot, qtot
end

function write_delta_readings(pf_sol,math)
    math_meas = deepcopy(math)

    for (b,bus) in pf_sol["bus"]
        bus["vmn2"] = []
        for t in setdiff(math["bus"][b]["terminals"],[4])
            push!(bus["vmn2"], (bus["vr"][t] - bus["vr"][4])^2 + (bus["vi"][t] - bus["vi"][4])^2)
        end
        bus["vmn"] = sqrt.(bus["vmn2"])
    math_meas["bus"][b]["terminals"] = setdiff(math["bus"][b]["terminals"], 4)
    end 

    for (l, load) in math["load"]
    if load["configuration"] == DELTA
        pf_sol["load"][l]["ptot"] = [sum(pf_sol["load"][l]["pd"])]
        pf_sol["load"][l]["qtot"] = [sum(pf_sol["load"][l]["qd"])]

        load_bus = pf_sol["bus"][string(load["load_bus"])]
        load_bus["vll"] = sqrt.([ (load_bus["vr"][x] - load_bus["vr"][y] )^2 + (load_bus["vi"][x] - load_bus["vi"][y] )^2 for (x,y) in [(1,2), (2,3), (3,1)]])
        end
    end
return math_meas
end




function kron_reduce_impedance(m::Matrix)
    # Check if the matrix is 4x4
    if size(m) != (4, 4)
        throw(DimensionMismatch("Input matrix must be 4x4."))
    end

    # Partition the matrix
    A = m[1:3, 1:3]  # 3x3 top-left block
    B = m[1:3, 4:4]  # 3x1 top-right block (column vector)
    C = m[4:4, 1:3]  # 1x3 bottom-left block (row vector)
    D_val = m[4, 4]    # 1x1 bottom-right block (scalar)

    # Check if D is invertible (non-zero for scalar)
    if isapprox(D_val, 0.0; atol=eps(real(ComplexF64)))
        error("Cannot perform Kron reduction: the pivot element m[4, 4] is zero or close to zero.")
        # Alternatively, depending on context, you might return an error or handle it differently.
        # For numerical stability, using inv(D) where D is a 1x1 matrix might be better
        # D = m[4:4, 4:4] # Represent D as a 1x1 matrix
        # if det(D) == 0 ... etc.
    end

    # Calculate the inverse of D (which is just 1/D_val for a scalar)
    D_inv_val = 1.0 / D_val

    # Calculate the Schur complement: A - B * D_inv * C
    # Note: B * D_inv_val * C performs matrix multiplication: (3x1) * scalar * (1x3) -> 3x3
    reduced_matrix = A - B * D_inv_val * C

    return reduced_matrix
end


"""
    get_sequence_components(m_phase::Matrix{ComplexF64})

Calculates the zero, positive, and negative sequence components of a 3x3
complex matrix representing phase quantities (e.g., impedance or admittance).

# Arguments
- `m_phase::Matrix{ComplexF64}`: The input 3x3 complex matrix in phase coordinates.

# Returns
- `Tuple{Matrix{ComplexF64}, ComplexF64, ComplexF64, ComplexF64}`: A tuple containing:
    1. The full 3x3 sequence matrix (M_seq).
    2. The zero sequence component (diagonal element M_seq[1, 1]).
    3. The positive sequence component (diagonal element M_seq[2, 2]).
    4. The negative sequence component (diagonal element M_seq[3, 3]).

# Throws
- `DimensionMismatch`: If the input matrix is not 3x3.


# Usage
Assume 'M_red' is the 3x3 matrix obtained from the previous Kron reduction step
Let's create a sample symmetrical 3x3 matrix for demonstration:
Z_phase = ComplexF64[
    10+5im   2+1im   2+1im;
    2+1im   10+5im  2+1im;
    2+1im   2+1im   10+5im
]
This represents a balanced system where off-diagonals are equal.

Or use a more general (unbalanced) example:
Z_phase = ComplexF64[
     10+5im   1+0.5im  2+1im;
     1+0.5im  12+6im   3+1.5im;
     2+1im    3+1.5im  11+5.5im
]

try
    # Calculate sequence components
    Z_seq_matrix, Z0, Z1, Z2 = calculate_sequence_components(Z_phase)

    println("Original Phase Matrix Z_phase:")
    display(Z_phase)

    println("\nSequence Matrix Z_seq:")
    display(round.(Z_seq_matrix; digits=4)) # Round for cleaner display

    println("\nSequence Components:")
    println("Zero Sequence (Z0): ", round(Z0; digits=4))
    println("Positive Sequence (Z1): ", round(Z1; digits=4))
    println("Negative Sequence (Z2): ", round(Z2; digits=4))

catch e
    println("Error: ", e)
end
"""
function get_sequence_components(m_phase::Matrix)
    # Check if the matrix is 3x3
    if size(m_phase) != (3, 3)
        throw(DimensionMismatch("Input matrix must be 3x3."))
    end

    # Define the complex operator 'a' (1 angle 120 degrees)
    α = exp(im * 2 * pi / 3) # cis(120°) or -0.5 + im*sqrt(3)/2

    # Define the Fortescue transformation matrix T
    T = ComplexF64[
        1  1      1;
        1  α^2    α;
        1  α      α^2
    ]

    # Define the inverse of the Fortescue transformation matrix T_inv
    # T_inv = (1/3) * T' conjugate transpose (Hermitian transpose)
    # Or explicitly:
    T_inv = (1/3) * ComplexF64[
        1  1      1;
        1  α      α^2;
        1  α^2    α
    ]

    # Calculate the sequence matrix: M_seq = T_inv * M_phase * T
    m_seq = T_inv * m_phase * T

    # Extract the diagonal components
    zero_sequence = m_seq[1, 1]
    positive_sequence = m_seq[2, 2]
    negative_sequence = m_seq[3, 3]

    return m_seq, zero_sequence, positive_sequence, negative_sequence
end




"""
    calculate_vuf(PF_RES, math)

Calculates the Voltage Unbalance Factor (VUF) for each bus in the power flow result.

VUF is defined as:
    VUF = |V₋ / V₊| × 100%
where:
- V₊ is the positive-sequence voltage
- V₋ is the negative-sequence voltage

It uses the phase-to-neutral complex voltages from the power flow solution, assumed
to be available at each bus under keys "1", "2", "3" (phases) and "4" (neutral)
after applying `dictify_solution!`.

Returns the modified PF_RES with VUF values added at each bus as:
    PF_RES["solution"]["bus"][bus_id]["vuf"] = <Float64>
"""
function calculate_vuf!(PF_RES, math)
    pf_sol = PF_RES["solution"]

    dictify_solution!(pf_sol, math)

    VUF_dict = Dict{String, Float64}()

    # get loaded buses 
    for (bus_id, bus_data) in pf_sol["bus"]
        voltages = bus_data["voltage"]

        # Extract phase-to-neutral voltages
        vag = haskey(voltages, "1") ? voltages["1"] : 0 + 0im
        vbg = haskey(voltages, "2") ? voltages["2"] : 0 + 0im
        vcg = haskey(voltages, "3") ? voltages["3"] : 0 + 0im
        vng = haskey(voltages, "4") ? voltages["4"] : 0 + 0im

        Va = vag - vng
        Vb = vbg - vng
        Vc = vcg - vng

        # Symmetrical component transformation
        V0, V1, V2 = calculate_symmetrical_components(Va, Vb, Vc)

        # Compute VUF
        vuf = abs(V2 / V1) * 100

        # Store VUF at the bus
        bus_data["vuf"] = vuf

        VUF_dict[bus_id] = vuf
    end

    return VUF_dict
end


function calculate_symmetrical_components(Va, Vb, Vc)
    # Operator 'α' = exp(j*2π/3)
    α = cis(2π/3)
    α² = α^2
    # Symmetrical component transformation
    V0 = (Va + Vb + Vc) / 3
    V1 = (Va + α  * Vb + α² * Vc) / 3  # Positive-sequence
    V2 = (Va + α² * Vb + α  * Vc) / 3  # Negative-sequence

    return V0, V1, V2
end


function calculate_symmetrical_components(bus_dict::Dict{String, Any})
    # Extract phase-to-neutral voltages
    vag = haskey(bus_dict["voltage"], "1") ? bus_dict["voltage"]["1"] : 0 + 0im
    vbg = haskey(bus_dict["voltage"], "2") ? bus_dict["voltage"]["2"] : 0 + 0im
    vcg = haskey(bus_dict["voltage"], "3") ? bus_dict["voltage"]["3"] : 0 + 0im
    vng = haskey(bus_dict["voltage"], "4") ? bus_dict["voltage"]["4"] : 0 + 0im

    Va = vag - vng
    Vb = vbg - vng
    Vc = vcg - vng

    return calculate_symmetrical_components(Va, Vb, Vc)
end

function plot_symmetrical_components(bus_dict::Dict{String, Any};     makie_backend=WGLMakie,
                    figure::Figure = nothing,
                    location::Tuple{Int, Int} = (1, 1),
                    fig_size=(800, 800),
                    kwargs...
                   )

    V0, V1, V2 = calculate_symmetrical_components(bus_dict)


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
                   title="Voltage Symemtrical Compoenents",
                   thetaticks=(radian_ticks, tick_labels),
                   kwargs...
                  )

    bus_phasor!(ax, eng, bus_id)

    return f, ax

end


function move_coords_eng_to_math!(eng, math)

    if !haskey(eng, "bus") || !haskey(math, "bus_lookup")
        error("The provided 'eng' or 'math' dictionary does not contain the required keys 'bus' or 'bus_lookup'.")
    end

    if !haskey(eng["bus"], "lat") || !haskey(eng["bus"], "lon")
        error("The 'eng' dictionary does not have 'lat' or 'lon' keys in the 'bus' sub-dictionary.")
    end
    for (b,bus) in eng["bus"]
        if b == "sourcebus"; continue end
        math["bus"]["$(math["bus_lookup"][b])"]["lat"] = bus["lat"]
        math["bus"]["$(math["bus_lookup"][b])"]["lon"] = bus["lon"]
    end

end

"""
    show_transformer_math_components(math; suppress_print::Bool=false)

Extract and display transformer mathematical components from a power system model.

# Arguments
- `math::Dict`: A dictionary containing the mathematical model representation with 
  "map" and element type keys (e.g., "bus", "branch", "transformer")
- `suppress_print::Bool=false`: If `true`, suppresses console output and only returns results

# Returns
- `Dict{String, Any}`: A nested dictionary organized by transformer name, containing:
  - `"buses"`: Dictionary of bus elements mapped to this transformer
  - `"branches"`: Dictionary of branch elements mapped to this transformer
  - `"transformers"`: Dictionary of transformer elements mapped to this transformer

# Description
This function identifies all transformers in the mathematical model that were converted 
from engineering model representation and extracts the associated mathematical components 
(buses, branches, transformers). Each transformer is identified by the presence of an 
"unmap_function" equal to "_map_math2eng_transformer!" in the mapping structure.

If `suppress_print=false`, displays formatted console output showing the transformation 
process and component details.

# Example
mapped_transformers = show_transformer_math_components(math, suppress_print=false)

"""
function show_transformer_math_components(math; suppress_print::Bool=false)
    results = Dict{String, Any}()   
    transformer_index_byname = findall(x->haskey(x, "unmap_function") && x["unmap_function"] == "_map_math2eng_transformer!", math["map"])
    for idx in transformer_index_byname
        trafo_name = math["map"][idx]["from"]
        
        if !suppress_print
            printstyled("👇 BEGIN Transformer ",color=:green)
            printstyled(trafo_name,color=:green, bold=true)
            printstyled(" 👇\n",color=:green)        
            display("in the MATHEMATICAL model transformer " * trafo_name*" was converted to: ")
        end
        
        results[trafo_name] = Dict{String, Any}("buses" => Dict{String, Any}(), "branches" => Dict{String, Any}(), "transformers" => Dict{String, Any}())
        
        for element in math["map"][idx]["to"]
            element_type = split(element, ".")[1]
            element_idx = split(element, ".")[2]
            
            if !suppress_print
                display("================================" * element_type * " "* element_idx *"================================")
                # display("Element Index: " * element_idx)    
                display(math[element_type][element_idx])
            end
            
            if element_type == "bus"
                results[trafo_name]["buses"][element_idx] = math[element_type][element_idx]
            elseif element_type == "branch"
                results[trafo_name]["branches"][element_idx] = math[element_type][element_idx]
            elseif element_type == "transformer"
                results[trafo_name]["transformers"][element_idx] = math[element_type][element_idx]
            end
        end
        
        if !suppress_print
            printstyled("👆 END Transformer ",color=:red)
            printstyled(trafo_name,color=:red, bold=true)
            printstyled(" 👆\n",color=:red)
        end
    end
    return results
end





"""
    add_degree_to_bus!(data)

Calculates the degree (number of connected branches) for each bus in the network data dictionary and adds it as a "degree" field to the respective bus dictionary.

# Arguments
- `data`: A dictionary containing PowerModels-style network data, specifically requiring "bus" and "branch" components.

# Effects
- Modifies `data["bus"]` in-place by adding a "degree" integer value to each bus entry.
"""
function add_degree_to_bus!(data)
    for (_, bus) in data["bus"]
        bus["degree"] = 0
        for (_, branch) in data["branch"]
            if branch["t_bus"] == bus["index"] || branch["f_bus"] == bus["index"]
                bus["degree"] += 1
            end
        end
    end
end




"""
    remove_all_superfluous_buses!(data::Dict)

Simplifies a PowerModelsDistribution MATHEMATICAL network dictionary by removing intermediate buses that have no load, no generation, and a degree of 2 or less (i.e., simple pass-through nodes).

# Arguments
- `data`: A dictionary containing MATHEMATICAL network data, specifically requiring "bus" and "branch" dictionaries.

# Functionality
The function performs the following steps:
1.  **Validation**: Checks that the dictionary is not a multi-network structure.
2.  **Bus Classification**: Identifies load buses and generation buses. It asserts that load buses have a degree of 1 (implying loads should not be on the main feeder line directly but connected via a spur).
3.  **Topology Analysis**: Calculates the degree of all buses.
4.  **Reduction Loop**:
    - Identifies candidate buses for deletion (non-load, non-gen, degree <= 2).
    - Iteratively merges branches connected to these buses.
    - Updates the resistance (`br_r`) and reactance (`br_x`) of the preserved branch by summing the values of the merged branches.
    - Updates connectivity (from/to bus indices) to bypass the deleted bus.
    - Deletes the superfluous bus and the redundant branch.
5.  **Reorientation**: Ensures that the branch connected to the reference bus (slack bus, type 3) is oriented such that `f_bus` is the reference bus.

# Returns
- The modified `data` dictionary with the topology simplified.
"""
function remove_all_superfluous_buses!(data::Dict)
    @assert !haskey(data, "nw") "Please use `remove_all_intermediate_buses_mn` for multinetwork data dicts like this one"
    load_buses = ["$(load["load_bus"])" for (_, load) in data["load"]]
    gen_buses = ["$(gen["gen_bus"])" for (_, gen) in data["gen"]]
    add_degree_to_bus!(data)
    for lb in load_buses @assert data["bus"][lb]["degree"] == 1 "Load $lb is on the main cable, add a small connection cable to that load" end
    to_delete = [b for (b, bus) in data["bus"] if (b ∉ union!(gen_buses, load_buses) && bus["degree"] <= 2)]
    for db in to_delete
        data["bus"][db]["adjacent_buses"] = []
        data["bus"][db]["inout_branches"] = []
        for (br, branch) in data["branch"]
            if branch["f_bus"] == parse(Int, db) || branch["t_bus"] == parse(Int, db) 
                push!(data["bus"][db]["inout_branches"], br)
                if branch["f_bus"] != parse(Int, db)
                    push!(data["bus"][db]["adjacent_buses"], "$(branch["f_bus"])")
                else
                    push!(data["bus"][db]["adjacent_buses"], "$(branch["t_bus"])")
                end
            end
        end
    end
    while !isempty(to_delete)
        for db in to_delete
            if any([b ∈ to_delete for b in data["bus"][db]["adjacent_buses"]])
                deletable_adj_bus = [b for b in data["bus"][db]["adjacent_buses"] if b ∈ to_delete][1]
                other_adj_bus = [b for b in data["bus"][db]["adjacent_buses"] if b != deletable_adj_bus][1]
                deletable_adj_bus_branches = data["bus"][deletable_adj_bus]["inout_branches"]
                delete_branch = first(intersect(Set(data["bus"][db]["inout_branches"]), Set(deletable_adj_bus_branches)))
                preserve_branch = [br for br in data["bus"][db]["inout_branches"] if br != delete_branch][1]
                
                Req = (data["branch"][preserve_branch]["br_r"] .+ data["branch"][delete_branch]["br_r"])
                Xeq = (data["branch"][preserve_branch]["br_x"] .+ data["branch"][delete_branch]["br_x"]) 
                data["branch"][preserve_branch]["br_r"] = Req
                data["branch"][preserve_branch]["br_x"] = Xeq
                
                data["branch"][preserve_branch]["f_bus"] = parse(Int64, other_adj_bus)
                data["branch"][preserve_branch]["t_bus"] = parse(Int64, deletable_adj_bus)
                data["bus"][deletable_adj_bus]["adjacent_buses"] = filter(x->x!=db, data["bus"][deletable_adj_bus]["adjacent_buses"])
                push!(data["bus"][deletable_adj_bus]["adjacent_buses"], other_adj_bus)
                data["bus"][deletable_adj_bus]["inout_branches"] = filter(x->x!=delete_branch, data["bus"][deletable_adj_bus]["inout_branches"])
                push!(data["bus"][deletable_adj_bus]["inout_branches"], preserve_branch)
                delete!(data["branch"], delete_branch)
            else
                delete_branch = data["bus"][db]["inout_branches"][1]
                preserve_branch = [br for br in data["bus"][db]["inout_branches"] if br != delete_branch][1]
                Req = (data["branch"][preserve_branch]["br_r"] .+ data["branch"][delete_branch]["br_r"])
                Xeq = (data["branch"][preserve_branch]["br_x"] .+ data["branch"][delete_branch]["br_x"]) 
                data["branch"][preserve_branch]["br_r"] = Req
                data["branch"][preserve_branch]["br_x"] = Xeq
                
                delete!(data["branch"], delete_branch)
                data["branch"][data["bus"][db]["inout_branches"][2]]["f_bus"] = parse(Int64, data["bus"][db]["adjacent_buses"][1])
                data["branch"][data["bus"][db]["inout_branches"][2]]["t_bus"] = parse(Int64, data["bus"][db]["adjacent_buses"][2])
            end
            delete!(data["bus"], db)
            to_delete = filter(x->x!=db, to_delete)
        end
    end
    # the lines below make sure that the orientation of the branch at the slack bus is from slack_bus to --> rest of feeder
    ref_bus = [bus["index"] for (_,bus) in data["bus"] if bus["bus_type"] == 3][1]
    ref_branch_fr = [b for (b, br) in data["branch"] if br["f_bus"] == ref_bus]
    if isempty(ref_branch_fr) 
        ref_branch_to = [b for (b, br) in data["branch"] if br["t_bus"] == ref_bus][1]
        f_bus = data["branch"][ref_branch_to]["f_bus"]
        data["branch"][ref_branch_to]["f_bus"] = ref_bus
        data["branch"][ref_branch_to]["t_bus"] = f_bus
    end
    return data
end


# function reduce_network_buses!(data::Dict)
#     @assert !haskey(data, "nw") "Multinetwork not supported in this function."

#     # 1. Pre-calculate bus connectivity (O(Branches))
#     # Maps bus_id -> list of branch_ids
#     bus_to_branches = Dict(b_idx => String[] for b_idx in keys(data["bus"]))
#     for (br_id, branch) in data["branch"]
#         push!(bus_to_branches["$(branch["f_bus"])"], br_id)
#         push!(bus_to_branches["$(branch["t_bus"])"], br_id)
#     end

#     # 2. Identify Load/Gen buses (Sets for O(1) lookup)
#     load_bus_ids = Set(["$(load["load_bus"])" for (_, load) in data["load"]])
#     gen_bus_ids = Set(["$(gen["gen_bus"])" for (_, gen) in data["gen"]])
    
#     # 3. Identify candidates for deletion (Degree == 2 and no Load/Gen)
#     to_delete = String[]
#     for (b_id, branch_list) in bus_to_branches
#         degree = length(branch_list)
        
#         # Validation: Load buses must be degree 1
#         if b_id in load_bus_ids
#             @assert degree == 1 "Load $b_id is on the main cable; degree is $degree"
#             continue
#         end

#         if degree == 2 && b_id ∉ gen_bus_ids
#             push!(to_delete, b_id)
#         end
#     end

#     # 4. Reduction Loop
#     for db in to_delete
#         # Get the two branches connected to this bus
#         br_ids = bus_to_branches[db]
#         if length(br_ids) != 2
#             continue # Already processed or modified by a neighbor's deletion
#         end

#         br1_id, br2_id = br_ids[1], br_ids[2]
#         br1 = data["branch"][br1_id]
#         br2 = data["branch"][br2_id]

#         # Find the two "far" buses (the ones not being deleted)
#         # We handle the fact that f_bus or t_bus could be the db
#         bus_a = "$(br1["f_bus"])" == db ? "$(br1["t_bus"])" : "$(br1["f_bus"])"
#         bus_b = "$(br2["f_bus"])" == db ? "$(br2["t_bus"])" : "$(br2["f_bus"])"

#         # Merge Impedance into br1 (Preserve br1, delete br2)
#         br1["br_r"] .+= br2["br_r"]
#         br1["br_x"] .+= br2["br_x"]

#         # Update br1 connectivity to bypass db
#         br1["f_bus"] = parse(Int, bus_a)
#         br1["t_bus"] = parse(Int, bus_b)

#         # Update the connectivity index for bus_b
#         # bus_b was connected to br2, now it's connected to br1
#         bus_to_branches[bus_b] = filter(id -> id != br2_id, bus_to_branches[bus_b])
#         push!(bus_to_branches[bus_b], br1_id)

#         # Remove the bus and the redundant branch
#         delete!(data["branch"], br2_id)
#         delete!(data["bus"], db)
#     end

#     # 5. Reorientation (Slack Bus)
#     ref_bus = [bus["index"] for (_, bus) in data["bus"] if bus["bus_type"] == 3][1]
#     for (br_id, br) in data["branch"]
#         if br["t_bus"] == ref_bus
#             # Flip orientation
#             f, t = br["f_bus"], br["t_bus"]
#             br["f_bus"], br["t_bus"] = t, f
#             break # Usually only one branch at the slack in radial feeders
#         end
#     end

#     return data
# end

"""
    reduce_network_buses!(data::Dict)

Simplifies a PowerModelsDistribution MATHEMATICAL network dictionary by removing intermediate buses that have no load, no generation, and a degree of 2 or less (i.e., simple pass-through nodes).

# Arguments
- `data`: A dictionary containing MATHEMATICAL network data, specifically requiring "bus" and "branch" dictionaries.

# Functionality
The function performs the following steps:
1.  **Validation**: Checks that the dictionary is not a multi-network structure.
2.  **Bus Classification**: Identifies load buses and generation buses. It asserts that load buses have a degree of 1 (implying loads should not be on the main feeder line directly but connected via a spur).
3.  **Topology Analysis**: Calculates the degree of all buses.
4.  **Reduction Loop**:
    - Identifies candidate buses for deletion (non-load, non-gen, degree <= 2).
    - Iteratively merges branches connected to these buses.
    - Updates the resistance (`br_r`) and reactance (`br_x`) of the preserved branch by summing the values of the merged branches.
    - Updates connectivity (from/to bus indices) to bypass the deleted bus.
    - Deletes the superfluous bus and the redundant branch.
5.  **Reorientation**: Ensures that the branch connected to the reference bus (slack bus, type 3) is oriented such that `f_bus` is the reference bus.

# Returns
- The modified `data` dictionary with the topology simplified.
"""
function reduce_network_buses!(data::Dict)
    @assert !haskey(data, "nw") "Multinetwork not supported in this function, apply before performing `make_multinetwork!`."

    # 1. Map Connectivity
    bus_to_branches = Dict(b => String[] for b in keys(data["bus"]))
    for (i, br) in data["branch"]
        push!(bus_to_branches["$(br["f_bus"])"], i)
        push!(bus_to_branches["$(br["t_bus"])"], i)
    end


    # Identify special_buses (to NOT delete it)
    special_buses = [bus["index"] for (_, bus) in data["bus"] if bus["bus_type"] != 1]
    @debug "sepcial_buses (not deletable): " * string(special_buses)

    # Identify buses connected to transformers (to NOT delete these)
    transformer_buses = Set{String}()
    if haskey(data, "transformer")
        for (_, tx) in data["transformer"]
            push!(transformer_buses, "$(tx["f_bus"])")
            push!(transformer_buses, "$(tx["t_bus"])")
        end
    end

    load_buses = Set(["$(l["load_bus"])" for (_, l) in data["load"]])
    gen_buses = Set(["$(g["gen_bus"])" for (_, g) in data["gen"]])

    # 2. Find Candidates (Degree 2, no load/gen, NOT a transformer terminal)
    to_delete = [b for (b, branches) in bus_to_branches if 
                 length(branches) == 2 &&  # if you want to remove hanging buses (loadless genless leaf buses) use <= 2
                 b ∉ load_buses && 
                 b ∉ gen_buses && 
                 b ∉ transformer_buses &&
                 b ∉ special_buses
                 ]

    # 3. Reduction Loop
    for db in to_delete
        br_ids = bus_to_branches[db]
        if length(br_ids) != 2 continue end # Re-check if neighbor already deleted

        id1, id2 = br_ids[1], br_ids[2]
        br1, br2 = data["branch"][id1], data["branch"][id2]

        # Determine connectivity: bus_a --(br1)-- db --(br2)-- bus_b
        bus_a = "$(br1["f_bus"])" == db ? "$(br1["t_bus"])" : "$(br1["f_bus"])"
        bus_b = "$(br2["f_bus"])" == db ? "$(br2["t_bus"])" : "$(br2["f_bus"])"

        # --- SHUNT TRANSFER LOGIC ---
        # If db had shunts (g_fr, b_fr, etc.) we must sum them into the preserved branch
        # We assume we are merging br2 into br1
        for key in ["g_fr", "g_to", "b_fr", "b_to"]
            # Sum shunts from the deleted bus 'db' into br1
            if "$(br2["f_bus"])" == db
                br1["g_to"] .+= br2["g_fr"]
                br1["b_to"] .+= br2["b_fr"]
            else
                br1["g_to"] .+= br2["g_to"]
                br1["b_to"] .+= br2["b_to"]
            end
        end

        # --- IMPEDANCE MERGE ---
        br1["br_r"] .+= br2["br_r"]
        br1["br_x"] .+= br2["br_x"]

        # Re-link br1 to bypass db
        br1["f_bus"] = parse(Int, bus_a)
        br1["t_bus"] = parse(Int, bus_b)

        # Update connectivity index for bus_b
        bus_to_branches[bus_b] = filter(id -> id != id2, bus_to_branches[bus_b])
        push!(bus_to_branches[bus_b], id1)

        delete!(data["branch"], id2)
        delete!(data["bus"], db)
    end

    # 5. Reorientation (Slack Bus)
    ref_bus = [bus["index"] for (_, bus) in data["bus"] if bus["bus_type"] == 3][1]
    data["settings"]
    for (br_id, br) in data["branch"]
        if br["t_bus"] == ref_bus
            # Flip orientation
            f, t = br["f_bus"], br["t_bus"]
            br["f_bus"], br["t_bus"] = t, f
            break # Usually only one branch at the slack in radial feeders
        end
    end
    return data
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


export bus_phasor, bus_phasor!, plot_bus_phasor, calculate_vuf!
export remove_all_superfluous_buses!, add_degree_to_bus!, reduce_network_buses!
end # module PMDUtils
