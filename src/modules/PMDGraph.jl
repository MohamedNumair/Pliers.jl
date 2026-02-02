"""
    PMDGraph

Internal sub-module providing plotting functions for power distribution network visualization.

This module re-exports functions from the main Pliers module for:
- Network tree visualization
- Coordinate-based network plotting
- Geographic map overlays
- Bus phasor diagram plotting

See the main Pliers module for function documentation.
"""
module PMDGraph

using ..Pliers
using ..PMDUtils
using ..PMDUtils: _is_eng


# plotting packages
using Makie
# using MakieCore
using CairoMakie
using WGLMakie
if Sys.iswindows()
    # using GLMakie
end

using GraphMakie
using GeoMakie
using Proj
using Graphs
using MetaGraphs





#=

░███    ░██ ░██████████ ░██████████░██       ░██   ░██████   ░█████████  ░██     ░██      ░██████  ░█████████     ░███    ░█████████  ░██     ░██ 
░████   ░██ ░██             ░██    ░██       ░██  ░██   ░██  ░██     ░██ ░██    ░██      ░██   ░██ ░██     ░██   ░██░██   ░██     ░██ ░██     ░██ 
░██░██  ░██ ░██             ░██    ░██  ░██  ░██ ░██     ░██ ░██     ░██ ░██   ░██      ░██        ░██     ░██  ░██  ░██  ░██     ░██ ░██     ░██ 
░██ ░██ ░██ ░█████████      ░██    ░██ ░████ ░██ ░██     ░██ ░█████████  ░███████       ░██  █████ ░█████████  ░█████████ ░█████████  ░██████████ 
░██  ░██░██ ░██             ░██    ░██░██ ░██░██ ░██     ░██ ░██   ░██   ░██   ░██      ░██     ██ ░██   ░██   ░██    ░██ ░██         ░██     ░██ 
░██   ░████ ░██             ░██    ░████   ░████  ░██   ░██  ░██    ░██  ░██    ░██      ░██  ░███ ░██    ░██  ░██    ░██ ░██         ░██     ░██ 
░██    ░███ ░██████████     ░██    ░███     ░███   ░██████   ░██     ░██ ░██     ░██      ░█████░█ ░██     ░██ ░██    ░██ ░██         ░██     ░██ 
                                                                                                                                                  
                                                                                                                                                  
                                                                                                                                                  
=#


# function _is_eng(graph::MetaDiGraph)

#     if haskey(first(graph.eprops).second, :branch_id)
#         return false
#     else
#         return true
#     end
#     error("I don't know if this graph is for a MATHEMATICAL or ENGINEERING model. I usually check if it has `:branch_id` or `:line_id` in the edge properties to tell which one it is.")
# end


"""
    create_network_graph(data::Dict{String,Any}, fallback_layout)

Create a MetaDiGraph representation of a power distribution network.

Automatically detects whether the input data is in ENGINEERING or MATHEMATICAL format
and calls the appropriate graph creation function.

# Arguments
- `data::Dict{String,Any}`: Network data dictionary (either ENGINEERING or MATHEMATICAL model).
- `fallback_layout`: Layout function to use if no coordinates are available.

# Returns
The result of either `create_network_graph_eng` or `create_network_graph_math` depending
on the data model type.

# See also
- [`create_network_graph_eng`](@ref)
- [`create_network_graph_math`](@ref)
"""
function create_network_graph(data::Dict{String,Any}, fallback_layout)
    _is_eng(data) ? create_network_graph_eng(data, fallback_layout) : create_network_graph_math(data, fallback_layout)
end

"""
    create_network_graph_eng(eng::Dict{String,Any}, fallback_layout) -> MetaDiGraph, Function, Dict{Symbol,Any}

Create a MetaDiGraph from an ENGINEERING model network.

Constructs a directed graph where vertices represent buses and edges represent
lines and transformers. Bus properties are enriched with connected loads, shunts,
and voltage sources.

# Arguments
- `eng::Dict{String,Any}`: Network data in ENGINEERING format.
- `fallback_layout`: Layout function to use if bus coordinates are not available.

# Returns
A tuple containing:
- `network_graph::MetaDiGraph`: The constructed network graph.
- `GraphLayout::Function`: Layout function for plotting (coordinate-based or fallback).
- `eng_sym::Dict{Symbol,Any}`: Copy of input data with keys converted to symbols.

# Examples
```julia
network_graph, layout, eng_sym = create_network_graph_eng(eng, GraphMakie.Spring())
```
"""
function create_network_graph_eng(eng::Dict{String,Any}, fallback_layout)
    eng_sym = convert_keys_to_symbols(deepcopy(eng))
    network_graph = MetaDiGraph()

    lons = []
    lats = []

    # Add bus keys as :bus_id and collect coordinates if present # and enrich the buses with loads and shunts connected to them
    for (b, bus) in eng_sym[:bus]
        bus[:bus_id] = b
        bus[:loads] = []
        bus[:shunts] = []
        # attach loads
        if haskey(eng_sym, :load)
            for (_, load) in eng_sym[:load]
                if Symbol(load[:bus]) == bus[:bus_id]
                    push!(bus[:loads], load)
                end
            end
        end
        # attach shunts
        if haskey(eng_sym, :shunt)
            for (_, shunt) in eng_sym[:shunt]
                if Symbol(shunt[:bus]) == bus[:bus_id]
                    push!(bus[:shunts], shunt)
                end
            end
        end

        # attach gens
        if haskey(eng_sym, :voltage_source)
            for (_, voltage_source) in eng_sym[:voltage_source]
                if Symbol(voltage_source[:bus]) == bus[:bus_id]
                    bus[:voltage_sources] = voltage_source
                end
            end
        end


        if haskey(bus, :lon) && haskey(bus, :lat)
            push!(lons, bus[:lon])
            push!(lats, bus[:lat])
        end
    end

    for (l, line) in eng_sym[:line]
        line[:line_id] = l
        
        haskey(line, :linecode) && (line[:linecodes] = eng_sym[:linecode][Symbol(line[:linecode])])

    end


    ### DETERMINE THE ROOT BUS
    sourcebus = eng_sym[:voltage_source][:source][:bus]

    # Determine source coordinates if available # and enrich the lines with the linecodes details
    lon_s, lat_s = nothing, nothing
    if length(lons) > 0
        source_line = findfirst(line -> line[:f_bus] == sourcebus, eng_sym[:line])
        if source_line !== nothing
            lon_s = eng_sym[:bus][Symbol(eng_sym[:line][source_line][:t_bus])][:lon]
            lat_s = eng_sym[:bus][Symbol(eng_sym[:line][source_line][:t_bus])][:lat]
        end
    end

    layouting_vector = []

    # Add `sourcebus` as the root
    if haskey(eng_sym[:bus], Symbol(sourcebus))
        add_vertex!(network_graph, eng_sym[:bus][Symbol(sourcebus)])
        if length(lons) > 0 && lon_s !== nothing && lat_s !== nothing
            push!(layouting_vector, (lon_s, lat_s))
        end
    else
        error("sourcebus not found in the bus data. Please add sourcebus to the bus data.")
    end

    # Add the rest of the buses
    for (_, bus) in eng_sym[:bus]
        if bus[:bus_id] != Symbol(sourcebus)
            add_vertex!(network_graph, bus)
            if haskey(bus, :lon) && haskey(bus, :lat)
                push!(layouting_vector, (bus[:lon], bus[:lat]))
            end
        end
    end

    # Use bus_id as the indexing property
    set_indexing_prop!(network_graph, :bus_id)

    # Add edges based on f_bus and t_bus
    for (_, line) in eng_sym[:line]
        f_bus = Symbol(line[:f_bus])
        t_bus = Symbol(line[:t_bus])
        f_vertex = network_graph[f_bus, :bus_id]
        t_vertex = network_graph[t_bus, :bus_id]
        add_edge!(network_graph, f_vertex, t_vertex, line)
    end

    # add transformers as edges based on f_bus and t_bus
    if haskey(eng_sym, :transformer)
        for (_, transformer) in eng_sym[:transformer]
            f_bus = Symbol(transformer[:bus][1])
            t_bus = Symbol(transformer[:bus][2])
            f_vertex = network_graph[f_bus, :bus_id]
            t_vertex = network_graph[t_bus, :bus_id]
            add_edge!(network_graph, f_vertex, t_vertex, transformer)
        end
    end
    # Decide on the layout
    if length(layouting_vector) > 1
        GraphLayout = _ -> layouting_vector
    else
        @warn "Note there were no coordinates found for plotting, the fallback (e.g. tree layout) layout will be used"
        GraphLayout = fallback_layout
    end

    return network_graph, GraphLayout, eng_sym
end


"""
    create_network_graph_math(math::Dict{String,Any}, fallback_layout) -> MetaDiGraph, Function, Dict{Symbol,Any}

Create a MetaDiGraph from a MATHEMATICAL model network.

Constructs a directed graph where vertices represent buses and edges represent
branches. Bus properties are enriched with connected loads, shunts, and generators.

# Arguments
- `math::Dict{String,Any}`: Network data in MATHEMATICAL format.
- `fallback_layout`: Layout function to use if bus coordinates are not available.

# Returns
A tuple containing:
- `network_graph::MetaDiGraph`: The constructed network graph.
- `GraphLayout::Function`: Layout function for plotting (coordinate-based or fallback).
- `math_sym::Dict{Symbol,Any}`: Copy of input data with keys converted to symbols.

# Examples
```julia
network_graph, layout, math_sym = create_network_graph_math(math, GraphMakie.Spring())
```
"""
function create_network_graph_math(math::Dict{String,Any}, fallback_layout)
    math_sym = convert_keys_to_symbols(deepcopy(math))
    network_graph = MetaDiGraph()

    lons = []
    lats = []

    # Add bus keys as :bus_id and collect coordinates if present # and enrich the buses with loads and shunts connected to them
    for (b, bus) in math_sym[:bus]
        bus[:bus_id] = b
        bus[:loads] = []
        bus[:shunts] = []
        if haskey(math_sym, :load)
            for (l, load) in math_sym[:load]
                if Symbol(load[:load_bus]) == bus[:bus_id]
                    push!(bus[:loads], load)
                end
            end
        end

        if haskey(math_sym, :shunt)
            for (_, shunt) in math_sym[:shunt]
                if Symbol(shunt[:shunt_bus]) == bus[:bus_id]
                    push!(bus[:shunts], shunt)
                end
            end
        end

        # gen 

        if haskey(math_sym, :gen)
            for (_, gen) in math_sym[:gen]
                if Symbol(gen[:gen_bus]) == bus[:bus_id]
                    bus[:gens] = gen
                end
            end
        end


        if haskey(bus, :lon) && haskey(bus, :lat)
            push!(lons, bus[:lon])
            push!(lats, bus[:lat])
        end
    end

    # Add branch keys as :branch_id 
    for (l, branch) in math_sym[:branch]
        branch[:branch_id] = l

    end



    # Determine source coordinates if available 
    lon_s, lat_s = nothing, nothing
    if length(lons) > 0

        sourcebus = findfirst(bus -> contains(bus["name"], "sourcebus"), math["bus"])
        display(sourcebus)
        virtual_branch = findfirst(branch -> contains(string(branch["f_bus"]), sourcebus), math["branch"])
        @warn "virtual branch"
        display(virtual_branch)
        display(math_sym[:branch][Symbol(virtual_branch)])

        gen_branch = findfirst(branch -> contains(string(branch["f_bus"]), sourcebus), math["branch"])
        @warn "gen branch"
        display(gen_branch)
        display(math_sym[:branch][Symbol(gen_branch)])

        #virtual_branch = findfirst(branch -> contains(branch["name"], "_virtual_branch.voltage_source.source"), math["branch"])


        if virtual_branch !== nothing
            # @warn "The virtual bus params"
            # display(math_sym[:branch][Symbol(virtual_branch)])  

            # @warn "The virtual branch params"
            # display(math_sym[:bus][Symbol(math_sym[:branch][Symbol(virtual_branch)][:t_bus])])


            lon_s = math_sym[:bus][Symbol(math_sym[:branch][Symbol(virtual_branch)][:t_bus])][:lon]
            lat_s = math_sym[:bus][Symbol(math_sym[:branch][Symbol(virtual_branch)][:t_bus])][:lat]
        end
    end

    layouting_vector = []

    # Add `sourcebus` as the root
    virtual_bus = findfirst(bus -> contains(bus[:name], "virtual_bus.voltage_source.source"), math_sym[:bus])

    if haskey(math_sym[:bus], virtual_bus)
        # Prefer the virtual bus own coordinates if they exist
        if haskey(math_sym[:bus][virtual_bus], :lon) && haskey(math_sym[:bus][virtual_bus], :lat)
            add_vertex!(network_graph, math_sym[:bus][virtual_bus])
            push!(layouting_vector, (math_sym[:bus][virtual_bus][:lon], math_sym[:bus][virtual_bus][:lat]))

        elseif lon_s !== nothing && lat_s !== nothing
            # fallback to inferred source branch target coords

            add_vertex!(network_graph, math_sym[:bus][virtual_bus])
            push!(layouting_vector, (lon_s, lat_s))

            add_vertex!(network_graph, math_sym[:bus][Symbol(math_sym[:branch][Symbol(virtual_branch)][:f_bus])])
            push!(layouting_vector, (lon_s, lat_s))


        end
    else
        error("sourcebus not found in the bus data. Please add sourcebus to the bus data.")
    end

    # Add the rest of the buses
    for (_, bus) in math_sym[:bus]
        if (bus[:bus_id] != virtual_bus) && (bus[:bus_id] != Symbol(math_sym[:branch][Symbol(virtual_branch)][:f_bus]))
            add_vertex!(network_graph, bus)
            if haskey(bus, :lon) && haskey(bus, :lat)
                push!(layouting_vector, (bus[:lon], bus[:lat]))
            end
        end
    end

    # Use bus_id as the indexing property
    set_indexing_prop!(network_graph, :bus_id)

    # Add edges based on f_bus and t_bus
    for (_, branch) in math_sym[:branch]
        f_bus = Symbol(branch[:f_bus])
        t_bus = Symbol(branch[:t_bus])
        f_vertex = network_graph[f_bus, :bus_id]
        t_vertex = network_graph[t_bus, :bus_id]
        add_edge!(network_graph, f_vertex, t_vertex, branch)
    end

    # Decide on the layout
    if length(layouting_vector) > 1
        # Ensure full coverage; if partial, fall back
        if length(layouting_vector) == nv(network_graph)
            GraphLayout = _ -> layouting_vector
        else
            @warn "Incomplete coordinate set ($(length(layouting_vector)) of $(nv(network_graph)) buses); using fallback layout"
            GraphLayout = fallback_layout
        end
    else
        @warn "Note there were no coordinates found for plotting, the fallback (e.g. tree layout) layout will be used"
        GraphLayout = fallback_layout
    end

    return network_graph, GraphLayout, math_sym
end


# fix MetaGraphs set_indexing_prop! function


"""
    set_indexing_prop!(g, prop)
    set_indexing_prop!(g, v, prop, val)

Make property `prop` into an indexing property. If any values for this property
are already set, each vertex must have unique values. Optionally, set the index
`val` for vertex `v`. Any vertices without values will be set to a default
("(prop)(v)").
"""
function set_indexing_prop!(g::AbstractMetaGraph, prop::Symbol; exclude=nothing)
    in(prop, g.indices) && return g.indices
    index_values = [g.vprops[v][prop] for v in keys(g.vprops) if haskey(g.vprops[v], prop)]
    length(index_values) != length(union(index_values)) && error("Cannot make $prop an index, duplicate values detected")
    index_values = Set(index_values)

    g.metaindex[prop] = Dict{Any,Integer}()
    for v in vertices(g)
        if !haskey(g.vprops, v) || !haskey(g.vprops[v], prop)
            val = default_index_value(v, prop, index_values, exclude=exclude)
            set_prop!(g, v, prop, val)
        end
        g.metaindex[prop][g.vprops[v][prop]] = v
    end
    push!(g.indices, prop)
    return g.indices
end


"""
    default_index_value(v, prop, index_values; exclude=nothing)

Provides a default index value for a vertex if no value currently exists. The default is a string: "\$prop\$i" where `prop` is the property name and `i` is the vertex number. If some other vertex already has this name, a randomized string is generated (though the way it is generated is deterministic).
"""
function default_index_value(v::Integer, prop::Symbol, index_values::Set{}; exclude=nothing)
    val = string(prop) * string(v)
    if in(val, index_values) || val == exclude
        seed!(v + hash(prop))
        val = randstring()
        @warn("'$(string(prop))$v' is already in index, setting ':$prop' for vertex $v to $val")
    end
    return val
end


"""
    get_graph_node(G, node)
    get_graph_node(G, node, key)

Get properties of a node (bus) from the network graph.

# Arguments
- `G::MetaDiGraph`: The network graph.
- `node`: Node identifier (bus ID).
- `key`: (Optional) Specific property key to retrieve.

# Returns
- Without `key`: Dictionary of all node properties.
- With `key`: Value of the specified property.

# Examples
```julia
props = get_graph_node(G, "bus1")
voltage = get_graph_node(G, "bus1", "vm")
```
"""
function get_graph_node(G, node)
    Gidx = G[Symbol(node), :bus_id]
    return props(G, Gidx)
end


function get_graph_node(G, node, key)
    Gidx = G[Symbol(node), :bus_id]
    return props(G, Gidx)[Symbol(key)]
end

"""
    get_graph_edge(G, edge_id)
    get_graph_edge(G, edge_id, key)

Get properties of an edge (line/branch) from the network graph.

# Arguments
- `G::MetaDiGraph`: The network graph.
- `edge_id`: Edge identifier (line ID for ENGINEERING, branch ID for MATHEMATICAL).
- `key`: (Optional) Specific property key to retrieve.

# Returns
- Without `key`: Dictionary of all edge properties.
- With `key`: Value of the specified property.

# Examples
```julia
props = get_graph_edge(G, "line1")
impedance = get_graph_edge(G, "line1", "rs")
```
"""
function get_graph_edge(G, edge_id)
    if _is_eng(G)
        Eid = findall(e -> e[:line_id] == Symbol(edge_id), G.eprops)
    else
        Eid = findall(e -> e[:branch_id] == Symbol(edge_id), G.eprops)
    end
    return G.eprops[Eid...]
end
function get_graph_edge(G, edge_id, key)
    if _is_eng(G)
        Eid = findall(e -> e[:line_id] == Symbol(edge_id), G.eprops)
    else
        Eid = findall(e -> e[:branch_id] == Symbol(edge_id), G.eprops)
    end
    return G.eprops[Eid...][Symbol(key)]
end

# edge_color = [get_graph_edge(G, edge)[:linecodes][:rs] for edge in edges(G)]


#=
░███    ░██ ░██████████ ░██████████░██       ░██   ░██████   ░█████████  ░██     ░██    ░█████████  ░██           ░██████   ░██████████
░████   ░██ ░██             ░██    ░██       ░██  ░██   ░██  ░██     ░██ ░██    ░██     ░██     ░██ ░██          ░██   ░██      ░██    
░██░██  ░██ ░██             ░██    ░██  ░██  ░██ ░██     ░██ ░██     ░██ ░██   ░██      ░██     ░██ ░██         ░██     ░██     ░██    
░██ ░██ ░██ ░█████████      ░██    ░██ ░████ ░██ ░██     ░██ ░█████████  ░███████       ░█████████  ░██         ░██     ░██     ░██    
░██  ░██░██ ░██             ░██    ░██░██ ░██░██ ░██     ░██ ░██   ░██   ░██   ░██      ░██         ░██         ░██     ░██     ░██    
░██   ░████ ░██             ░██    ░████   ░████  ░██   ░██  ░██    ░██  ░██    ░██     ░██         ░██          ░██   ░██      ░██    
░██    ░███ ░██████████     ░██    ░███     ░███   ░██████   ░██     ░██ ░██     ░██    ░██         ░██████████   ░██████       ░██    
                                                                                                                                                                                                                                                                                                                                                                                                                     
=#



#=
 _______   __              __      __      __                     
/       \ /  |            /  |    /  |    /  |                    
$$$$$$$  |$$ |  ______   _$$ |_  _$$ |_   $$/  _______    ______  
$$ |__$$ |$$ | /      \ / $$   |/ $$   |  /  |/       \  /      \ 
$$    $$/ $$ |/$$$$$$  |$$$$$$/ $$$$$$/   $$ |$$$$$$$  |/$$$$$$  |
$$$$$$$/  $$ |$$ |  $$ |  $$ | __ $$ | __ $$ |$$ |  $$ |$$ |  $$ |
$$ |      $$ |$$ \__$$ |  $$ |/  |$$ |/  |$$ |$$ |  $$ |$$ \__$$ |
$$ |      $$ |$$    $$/   $$  $$/ $$  $$/ $$ |$$ |  $$ |$$    $$ |
$$/       $$/  $$$$$$/     $$$$/   $$$$/  $$/ $$/   $$/  $$$$$$$ |
                                                        /  \__$$ |
                                                        $$    $$/ 
                                                         $$$$$$/  
=#



#=
 _______   __              __            ________                                __      __                     
/       \ /  |            /  |          /        |                              /  |    /  |                    
$$$$$$$  |$$ |  ______   _$$ |_         $$$$$$$$/__    __  _______    _______  _$$ |_   $$/   ______   _______  
$$ |__$$ |$$ | /      \ / $$   |        $$ |__  /  |  /  |/       \  /       |/ $$   |  /  | /      \ /       \ 
$$    $$/ $$ |/$$$$$$  |$$$$$$/         $$    | $$ |  $$ |$$$$$$$  |/$$$$$$$/ $$$$$$/   $$ |/$$$$$$  |$$$$$$$  |
$$$$$$$/  $$ |$$ |  $$ |  $$ | __       $$$$$/  $$ |  $$ |$$ |  $$ |$$ |        $$ | __ $$ |$$ |  $$ |$$ |  $$ |
$$ |      $$ |$$ \__$$ |  $$ |/  |      $$ |    $$ \__$$ |$$ |  $$ |$$ \_____   $$ |/  |$$ |$$ \__$$ |$$ |  $$ |
$$ |      $$ |$$    $$/   $$  $$/       $$ |    $$    $$/ $$ |  $$ |$$       |  $$  $$/ $$ |$$    $$/ $$ |  $$ |
$$/       $$/  $$$$$$/     $$$$/        $$/      $$$$$$/  $$/   $$/  $$$$$$$/    $$$$/  $$/  $$$$$$/  $$/   $$/ 



=#


# Coord Transform

transSpanishTo4326 = Proj.Transformation("EPSG:3042", "EPSG:4326", always_xy=true) # ETRS89 / UTM zone 30N to WGS84
#transSpanishTo4326 = Proj.Transformation("EPSG:25830", "EPSG:4326", always_xy=true) # ETRS89 / UTM zone 30N to WGS84
transBritishto4326 = Proj.Transformation("EPSG:27700", "EPSG:4326", always_xy=true) # British National Grid OSGB36 to WGS84



"""
    network_graph_plot(network_graph::MetaDiGraph; layout=GraphMakie.Spring(), figure_size=(1000, 1200), node_color=:black, node_size=automatic, node_marker=automatic, node_strokewidth=automatic, show_node_labels=false, show_edge_labels=false, edge_color=:black, elabels_color=:black, elabels_fontsize=10, tangents=((0,-1),(0,-1)), arrow_show=false, arrow_marker='➤', arrow_size=50, arrow_shift=0.5)

Plots the given network graph using the specified layout and styling options.

# Arguments
- `network_graph::MetaDiGraph`: The network graph to be plotted.

# Keyword Arguments
- `layout`: The layout algorithm to use for positioning the nodes (default: `GraphMakie.Spring()`).
- `figure_size`: The size of the figure in pixels (default: `(1000, 1200)`).
- `node_color`: The color of the nodes (default: `:black`).
- `node_size`: The size of the nodes (default: `automatic`).
- `node_marker`: The marker style for the nodes (default: `automatic`).
- `node_strokewidth`: The stroke width of the nodes (default: `automatic`).
- `show_node_labels`: Whether to show labels on the nodes (default: `false`).
- `show_edge_labels`: Whether to show labels on the edges (default: `false`).
- `edge_color`: The color of the edges (default: `:black`).
- `elabels_color`: The color of the edge labels (default: `:black`).
- `elabels_fontsize`: The font size of the edge labels (default: `10`).
- `tangents`: The tangents for the edges (default: `((0,-1),(0,-1))`).
- `arrow_show`: Whether to show arrows on the edges (default: `false`).
- `arrow_marker`: The marker style for the arrows (default: `'➤'`).
- `arrow_size`: The size of the arrows (default: `50`).
- `arrow_shift`: The shift of the arrows along the edges (default: `0.5`).

# Returns
- A plot object representing the network graph.

# Example
```julia
f, ax, p = network_graph_plot(network_graph; layout=GraphMakie.Spring(), figure_size=(1000, 1200), node_color=:black, node_size=automatic, node_marker=automatic, node_strokewidth=automatic, show_node_labels=false, show_edge_labels=false, edge_color=:black, elabels_color=:black, elabels_fontsize=10, tangents=((0,-1),(0,-1)), arrow_show=false, arrow_marker='➤', arrow_size=50, arrow_shift=0.5)
```
"""
function network_graph_plot(
    network_graph::MetaDiGraph;

    # figure
    makie_backend=WGLMakie,
    layout=GraphMakie.Spring(), #layout=GraphMakie.Spring(),
    figure_size=(1000, 1200), #resolution

    #nodes
    nlabels=nothing,
    ilabels=nothing,
    node_color=:black,
    node_size=automatic,
    node_marker=automatic,
    node_strokewidth=automatic,
    show_node_labels=false, #nlablels
    nlabels_fontsize=12,
    #edges
    elabels=nothing,
    show_edge_labels=false, #elabels
    edge_color=:black,
    elabels_color=:black,
    elabels_fontsize=12,
    tangents=nothing,

    # arrow
    arrow_show=false,
    arrow_marker='➤',
    arrow_size=12,
    arrow_shift=0.5,
    kwargs...)

    makie_backend.activate!()
    #labels_theme = default_theme(makie_backend.Scene(), Makie.Text)
    #elabels_fontsize = isnothing(elabels_fontsize) ? labels_theme.fontsize : elabels_fontsize
    @debug "The used layout is: $layout"
    return graphplot(
        network_graph;
        layout=layout,
        resolution=figure_size,

        #nodes
        nlabels=nlabels,
        ilabels=ilabels,
        node_color=node_color,
        node_size=node_size,
        node_marker=node_marker,
        node_strokewidth=node_strokewidth,
        show_node_labels=show_node_labels,
        nlabels_fontsize=nlabels_fontsize,

        #edges
        elabels=elabels,
        show_edge_labels=show_edge_labels,
        edge_color=edge_color,
        elabels_color=elabels_color,
        elabels_fontsize=elabels_fontsize,
        tangents=tangents,

        # arrow
        arrow_show=arrow_show,
        arrow_marker=arrow_marker,
        arrow_size=arrow_size,
        arrow_shift=arrow_shift,

        # edges are linesegments
        edge_plottype=:linesegments, kwargs...
    )
end


function network_graph_plot!(
    network_graph::MetaDiGraph;

    # figure
    makie_backend=WGLMakie,
    layout=GraphMakie.Spring(),
    figure_size=(1000, 1200), #resolution

    #nodes
    nlabels=nothing,
    ilabels=nothing,
    node_color=:black,
    node_size=automatic,
    node_marker=automatic,
    node_strokewidth=automatic,
    show_node_labels=false, #nlablels

    #edges
    elabels=nothing,
    show_edge_labels=false, #elabels
    edge_color=:black,
    elabels_color=:black,
    #elabels_fontsize =  nothing,
    #tangents =((0,-1),(0,-1)),
    tangents=nothing,

    # arrow
    arrow_show=false,
    arrow_marker='➤',
    arrow_size=12,
    arrow_shift=0.5,
    kwargs...)

    #makie_backend.activate!()
    #labels_theme = default_theme(makie_backend.Scene(), Makie.Text)
    #elabels_fontsize = isnothing(elabels_fontsize) ? labels_theme.fontsize : elabels_fontsize

    add_graph = graphplot!(
        network_graph;
        layout=layout,
        resolution=figure_size,

        #nodes
        nlabels=nlabels,
        ilabels=ilabels,
        node_color=node_color,
        node_size=node_size,
        node_marker=node_marker,
        node_strokewidth=node_strokewidth,
        show_node_labels=show_node_labels,

        #edges
        elabels=elabels,
        show_edge_labels=show_edge_labels,
        edge_color=edge_color,
        elabels_color=elabels_color,
        #elabels_fontsize=elabels_fontsize,
        tangents=tangents,

        # arrow
        arrow_show=arrow_show,
        arrow_marker=arrow_marker,
        arrow_size=arrow_size,
        arrow_shift=arrow_shift,


        # edges are linesegments
        edge_plottype=:linesegments,
        kwargs...
    )
    translate!(add_graph, 0, 0, 100)
end



function network_graph_plot!(
    ax::Axis,
    network_graph::MetaDiGraph;

    # figure
    makie_backend=WGLMakie,
    layout=GraphMakie.Spring(),
    figure_size=(1000, 1200), #resolution

    #nodes
    nlabels=nothing,
    ilabels=nothing,
    node_color=:black,
    node_size=automatic,
    node_marker=automatic,
    node_strokewidth=automatic,
    show_node_labels=false, #nlablels

    #edges
    elabels=nothing,
    show_edge_labels=false, #elabels
    edge_color=:black,
    elabels_color=:black,
    #elabels_fontsize =  nothing,
    tangents=nothing,

    # arrow
    arrow_show=false,
    arrow_marker='➤',
    arrow_size=12,
    arrow_shift=0.5,
    kwargs...)

    #makie_backend.activate!()
    #labels_theme = default_theme(makie_backend.Scene(), Makie.Text)
    #elabels_fontsize = isnothing(elabels_fontsize) ? labels_theme.fontsize : elabels_fontsize

    add_graph = graphplot!(
        ax,
        network_graph;
        layout=layout,
        resolution=figure_size,

        #nodes
        nlabels=nlabels,
        ilabels=ilabels,
        node_color=node_color,
        node_size=node_size,
        node_marker=node_marker,
        node_strokewidth=node_strokewidth,
        show_node_labels=show_node_labels,

        #edges
        elabels=elabels,
        show_edge_labels=show_edge_labels,
        edge_color=edge_color,
        elabels_color=elabels_color,
        #elabels_fontsize=elabels_fontsize,
        tangents=tangents,

        # arrow
        arrow_show=arrow_show,
        arrow_marker=arrow_marker,
        arrow_size=arrow_size,
        arrow_shift=arrow_shift,


        # edges are linesegments
        edge_plottype=:linesegments,
        kwargs...
    )
    translate!(add_graph, 0, 0, 100)
end


"""
    network_graph_map_plot(
        network_graph::MetaDiGraph,
        GraphLayout::Function;
        tiles_provider=TileProviders.Google(:satelite),
        zoom_lon=0.0942,
        zoom_lat=0.0942,
        makie_backend=WGLMakie,
        figure_size=(1000, 1200),
        nlabels=nothing,
        ilabels=nothing,
        node_color=:black,
        node_size=automatic,
        node_marker=automatic,
        node_strokewidth=automatic,
        show_node_labels=false,
        elabels=nothing,
        show_edge_labels=false,
        edge_color=:black,
        elabels_color=:black,
        elabels_fontsize=10,
        tangents=((0,-1),(0,-1)),
        arrow_show=false,
        arrow_marker='➤',
        arrow_size=12,
        arrow_shift=0.5,
        kwargs...
    )

Plots a network graph on a map using the specified layout and visual properties.

# Arguments
- `network_graph::MetaDiGraph`: The network graph to be plotted.
- `GraphLayout::Function`: The layout function for positioning the nodes.
- `tiles_provider`: The tile provider for the map background. Default is `TileProviders.Google(:satelite)`.
- `zoom_lon`: The longitudinal zoom level. Default is `0.0942`.
- `zoom_lat`: The latitudinal zoom level. Default is `0.0942`.
- `makie_backend`: The Makie backend to use for plotting. Default is `WGLMakie`.
- `figure_size`: The size of the figure in pixels. Default is `(1000, 1200)`.
- `nlabels`: Node labels. Default is `nothing`.
- `ilabels`: Internal labels. Default is `nothing`.
- `node_color`: Color of the nodes. Default is `:black`.
- `node_size`: Size of the nodes. Default is `automatic`.
- `node_marker`: Marker style for the nodes. Default is `automatic`.
- `node_strokewidth`: Stroke width of the nodes. Default is `automatic`.
- `show_node_labels`: Whether to show node labels. Default is `false`.
- `elabels`: Edge labels. Default is `nothing`.
- `show_edge_labels`: Whether to show edge labels. Default is `false`.
- `edge_color`: Color of the edges. Default is `:black`.
- `elabels_color`: Color of the edge labels. Default is `:black`.
- `elabels_fontsize`: Font size of the edge labels. Default is `10`.
- `tangents`: Tangents for the edges. Default is `((0,-1),(0,-1))`.
- `arrow_show`: Whether to show arrows on the edges. Default is `false`.
- `arrow_marker`: Marker style for the arrows. Default is `'➤'`.
- `arrow_size`: Size of the arrows. Default is `12`.
- `arrow_shift`: Shift of the arrows. Default is `0.5`.
- `kwargs`: Additional keyword arguments.

# Returns
- A map object with the plotted network graph.
"""
function network_graph_map_plot(
    network_graph::MetaDiGraph,
    GraphLayout::Function;
    # map 
    tiles_provider=TileProviders.Google(:satelite), # :roadmap, :satelite, :terrain, :hybrid
    zoom_lon=0.0942,
    zoom_lat=0.0942,
    # figure
    makie_backend=WGLMakie,
    figure_size=(1000, 1200), #resolution

    #nodes
    nlabels=nothing,
    ilabels=nothing,
    node_color=:black,
    node_size=automatic,
    node_marker=automatic,
    node_strokewidth=automatic,
    show_node_labels=false, #nlablels

    #edges
    elabels=nothing,
    show_edge_labels=false, #elabels
    edge_color=:black,
    elabels_color=:black,
    elabels_fontsize=10,
    tangents=nothing,

    # arrow
    arrow_show=false,
    arrow_marker='➤',
    arrow_size=12,
    arrow_shift=0.5,
    kwargs...)

    map_layout = GraphLayout(1)
    center_lon = mean(lonslats[1] for lonslats in map_layout)
    center_lat = mean(lonslats[2] for lonslats in map_layout)
    max_lon = maximum(lonslats[1] for lonslats in map_layout)
    max_lat = minimum(lonslats[2] for lonslats in map_layout)
    cent_to_max_lon = abs(max_lon - center_lon) * 2
    cent_to_max_lat = abs(max_lat - center_lat) * 2

    zoom_lon == 0.0942 ? zoom_lon = cent_to_max_lon : zoom_lon = zoom_lon
    zoom_lat == 0.0942 ? zoom_lat = cent_to_max_lat : zoom_lat = zoom_lat

    makie_backend.activate!()
    map_window_coords = Rect2f(center_lon - zoom_lon / 2, center_lat - zoom_lat / 2, zoom_lon, zoom_lat)


    # m = Tyler.Map(map_window_coords; provider=tiles_provider, crs=Tyler.wgs84)
    # hidedecorations!(m.axis)
    # hidespines!(m.axis)





    graph = graphplot!(
        network_graph;
        layout=GraphLayout,
        resolution=figure_size,

        #nodes
        nlabels=nlabels,
        ilabels=ilabels,
        node_color=node_color,
        node_size=node_size,
        node_marker=node_marker,
        node_strokewidth=node_strokewidth,
        show_node_labels=show_node_labels,

        #edges
        elabels=elabels,
        show_edge_labels=show_edge_labels,
        edge_color=edge_color,
        elabels_color=elabels_color,
        elabels_fontsize=elabels_fontsize,
        tangents=tangents,

        # arrow
        arrow_show=arrow_show,
        arrow_marker=arrow_marker,
        arrow_size=arrow_size,
        arrow_shift=arrow_shift,


        # edges are linesegments
        edge_plottype=:linesegments,
        kwargs...
    )

    translate!(graph, 0, 0, 100)

    return m #, m.axis, graph
end

#=
 _______   __              __            _______                       __                               
/       \ /  |            /  |          /       \                     /  |                              
$$$$$$$  |$$ |  ______   _$$ |_         $$$$$$$  |  ______    _______ $$/   ______    ______    _______ 
$$ |__$$ |$$ | /      \ / $$   |        $$ |__$$ | /      \  /       |/  | /      \  /      \  /       |
$$    $$/ $$ |/$$$$$$  |$$$$$$/         $$    $$< /$$$$$$  |/$$$$$$$/ $$ |/$$$$$$  |/$$$$$$  |/$$$$$$$/ 
$$$$$$$/  $$ |$$ |  $$ |  $$ | __       $$$$$$$  |$$    $$ |$$ |      $$ |$$ |  $$ |$$    $$ |$$      \ 
$$ |      $$ |$$ \__$$ |  $$ |/  |      $$ |  $$ |$$$$$$$$/ $$ \_____ $$ |$$ |__$$ |$$$$$$$$/  $$$$$$  |
$$ |      $$ |$$    $$/   $$  $$/       $$ |  $$ |$$       |$$       |$$ |$$    $$/ $$       |/     $$/ 
$$/       $$/  $$$$$$/     $$$$/        $$/   $$/  $$$$$$$/  $$$$$$$/ $$/ $$$$$$$/   $$$$$$$/ $$$$$$$/  
                                                                          $$ |                          
                                                                          $$ |                          
                                                                          $$/                           
=#

"""
    plot_network_tree(eng::Dict{String,Any}; makie_backend=WGLMakie)

Plots a network tree based on the given engineering model dictionary `eng`.

# Arguments
- `eng::Dict{String,Any}`: A dictionary containing the engineering model data.
- `makie_backend`: The backend to use for plotting. Defaults to `WGLMakie`.
- `network_graph::MetaDiGraph`: The network graph to be plotted.
- `GraphLayout::Function`: The layout function for positioning the nodes.
- `tiles_provider`: The tile provider for the map background. Default is `TileProviders.Google(:satelite)`.
- `zoom_lon`: The longitudinal zoom level. Default is `0.0942`.
- `zoom_lat`: The latitudinal zoom level. Default is `0.0942`.
- `makie_backend`: The Makie backend to use for plotting. Default is `WGLMakie`.
- `figure_size`: The size of the figure in pixels. Default is `(1000, 1200)`.
- `nlabels`: Node labels. Default is `nothing`.
- `ilabels`: Internal labels. Default is `nothing`.
- `node_color`: Color of the nodes. Default is `:black`.
- `node_size`: Size of the nodes. Default is `automatic`.
- `node_marker`: Marker style for the nodes. Default is `automatic`.
- `node_strokewidth`: Stroke width of the nodes. Default is `automatic`.
- `show_node_labels`: Whether to show node labels. Default is `false`.
- `elabels`: Edge labels. Default is `nothing`.
- `show_edge_labels`: Whether to show edge labels. Default is `false`.
- `edge_color`: Color of the edges. Default is `:black`.
- `elabels_color`: Color of the edge labels. Default is `:black`.
- `elabels_fontsize`: Font size of the edge labels. Default is `10`.
- `tangents`: Tangents for the edges. Default is `((0,-1),(0,-1))`.
- `arrow_show`: Whether to show arrows on the edges. Default is `false`.
- `arrow_marker`: Marker style for the arrows. Default is `'➤'`.
- `arrow_size`: Size of the arrows. Default is `12`.
- `arrow_shift`: Shift of the arrows. Default is `0.5`.
- `kwargs`: Additional keyword arguments.

# Returns
- A tuple `(f, ax, p)` where `f` is the figure, `ax` is the axis, and `p` is the plot object.

# Details
1. Activates the specified Makie backend.
2. Converts the keys of the engineering model dictionary to symbols.
3. Creates a `MetaDiGraph` to represent the network graph.
4. Adds bus keys as `:bus_id` and line keys as `:line_id`.
5. Adds the `sourcebus` as the root vertex of the network graph.
6. Adds the rest of the buses to the network graph.
7. Sets the indexing property of the network graph to `:bus_id`.
8. Adds edges to the network graph based on the `f_bus` and `t_bus` indices.
9. Plots the network graph using `graphplot` with labels for nodes and edges.
10. Hides decorations and spines of the plot axis.

# Errors
- Throws an error if `sourcebus` is not found in the bus data.
"""
function plot_network_tree(
    data::Dict{String,Any};
    figure_size=(1000, 1200),
    show_node_labels=false,
    show_edge_labels=false,
    layout=GraphMakie.Spring(),
    edge_labels_type=:line_id,
    phase="1",
    kwargs...
)
    # Create the network meta graph 
    network_graph, _, _ = create_network_graph(data, layout)

    # Handle labels if required
    nlabels = show_node_labels ? _write_nlabels(network_graph, data) : nothing
    elabels = show_edge_labels ? edge_labels_type == :line_id ? _write_line_id_elabels(network_graph, data) : _write_results_elabels(network_graph, data, phase) : nothing

    # FORCED NODE FORMATTING:
    # 1. Identify the sourcebus
    if _is_eng(data)
        _decorate_nodes!(network_graph, data)
        node_color = [props(network_graph, i)[:node_color] for i in 1:nv(network_graph)]
        node_marker = [props(network_graph, i)[:node_marker] for i in 1:nv(network_graph)]
        node_size = [props(network_graph, i)[:marker_size] for i in 1:nv(network_graph)]
        _decorate_edges(network_graph, data)
        edge_color = [get_prop(network_graph, e, :edge_color) for e in edges(network_graph)]
        arrow_show = [get_prop(network_graph, e, :arrow_show) for e in edges(network_graph)]
        arrow_marker = [get_prop(network_graph, e, :arrow_marker) for e in edges(network_graph)]
        arrow_size = [get_prop(network_graph, e, :arrow_size) for e in edges(network_graph)]
        arrow_shift = [get_prop(network_graph, e, :arrow_shift) for e in edges(network_graph)]
    else
        node_color = :black
        node_marker = :circle
        node_size = 10
        edge_color = :black
        arrow_show = false
        arrow_marker = '➤'
        arrow_size = 12
        arrow_shift = 0.5
    end
    # plot and return the network 
    return network_graph_plot(
        network_graph;
        layout=layout,
        figure_size=figure_size,
        show_node_labels=show_node_labels,
        nlabels=nlabels,
        show_edge_labels=show_edge_labels,
        elabels=elabels,
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
end






"""
    plot_network_tree(dss::String; kwargs...)

Plots a network tree from a given DSS file.

# Arguments
- `dss::String`: The path to the DSS file containing the network data.
- `makie_backend`: The Makie backend to use for plotting. Defaults to `WGLMakie`.

# Description
This function parses the DSS file to create an engineering model of the network, transforms loops in the model, and converts keys to symbols. It then constructs a network graph using `MetaDiGraph`, adds vertices for each bus, and sets the `sourcebus` as the root. Edges are added based on the `f_bus` and `t_bus` indices of the lines. Finally, it plots the network graph using `graphplot` with a specified layout and labels for nodes and edges.

# Returns
A plot of the network graph.
"""
function plot_network_tree(dss::String; kwargs...)
    eng = PowerModelsDistribution.parse_file(dss)
    plot_network_tree(eng, kwargs...)
end

function plot_network_tree!(
    ax::Axis,
    data::Dict{String,Any};
    figure_size=(1000, 1200),
    show_node_labels=false,
    show_edge_labels=false,
    layout=GraphMakie.Spring(),
    edge_labels_type=:line_id,
    phase="1",
    kwargs...
)
    # Create the network meta graph 
    network_graph, _, _ = create_network_graph(data, layout)

    # Handle labels if required
    nlabels = show_node_labels ? _write_nlabels(network_graph, data) : nothing
    elabels = show_edge_labels ? edge_labels_type == :line_id ? _write_line_id_elabels(network_graph, data) : _write_results_elabels(network_graph, data, phase) : nothing

    # FORCED NODE FORMATTING:
    # 1. Identify the sourcebus
    if _is_eng(data)
        _decorate_nodes!(network_graph, data)
        node_color = [props(network_graph, i)[:node_color] for i in 1:nv(network_graph)]
        node_marker = [props(network_graph, i)[:node_marker] for i in 1:nv(network_graph)]
        node_size = [props(network_graph, i)[:marker_size] for i in 1:nv(network_graph)]
        _decorate_edges(network_graph, data)
        edge_color = [get_prop(network_graph, e, :edge_color) for e in edges(network_graph)]
        arrow_show = [get_prop(network_graph, e, :arrow_show) for e in edges(network_graph)]
        arrow_marker = [get_prop(network_graph, e, :arrow_marker) for e in edges(network_graph)]
        arrow_size = [get_prop(network_graph, e, :arrow_size) for e in edges(network_graph)]
        arrow_shift = [get_prop(network_graph, e, :arrow_shift) for e in edges(network_graph)]
    else
        node_color = :black
        node_marker = :circle
        node_size = 10
        edge_color = :black
        arrow_show = false
        arrow_marker = '➤'
        arrow_size = 12
        arrow_shift = 0.5
    end
    # plot and return the network 
    return network_graph_plot!(
        ax,
        network_graph;
        layout=layout,
        figure_size=figure_size,
        show_node_labels=show_node_labels,
        nlabels=nlabels,
        show_edge_labels=show_edge_labels,
        elabels=elabels,
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
end


"""
    plot_network_coords(eng::Dict{String,Any}; show_node_labels=false, show_edge_labels=false, fallback_layout=GraphMakie.Spring(), kwargs...)

Plot a network graph with optional node and edge labels.

# Arguments
- `eng::Dict{String,Any}`: The engineering data used to create the network graph.
- `show_node_labels::Bool`: Whether to display labels for the nodes. Default is `false`.
- `show_edge_labels::Bool`: Whether to display labels for the edges. Default is `false`.
- `fallback_layout`: The layout algorithm to use if no specific layout is provided. Default is `GraphMakie.Spring()`.
- `network_graph::MetaDiGraph`: The network graph to be plotted.
- `GraphLayout::Function`: The layout function for positioning the nodes.
- `tiles_provider`: The tile provider for the map background. Default is `TileProviders.Google(:satelite)`.
- `zoom_lon`: The longitudinal zoom level. Default is `0.0942`.
- `zoom_lat`: The latitudinal zoom level. Default is `0.0942`.
- `makie_backend`: The Makie backend to use for plotting. Default is `WGLMakie`.
- `figure_size`: The size of the figure in pixels. Default is `(1000, 1200)`.
- `nlabels`: Node labels. Default is `nothing`.
- `ilabels`: Internal labels. Default is `nothing`.
- `node_color`: Color of the nodes. Default is `:black`.
- `node_size`: Size of the nodes. Default is `automatic`.
- `node_marker`: Marker style for the nodes. Default is `automatic`.
- `node_strokewidth`: Stroke width of the nodes. Default is `automatic`.
- `show_node_labels`: Whether to show node labels. Default is `false`.
- `elabels`: Edge labels. Default is `nothing`.
- `show_edge_labels`: Whether to show edge labels. Default is `false`.
- `edge_color`: Color of the edges. Default is `:black`.
- `elabels_color`: Color of the edge labels. Default is `:black`.
- `elabels_fontsize`: Font size of the edge labels. Default is `10`.
- `tangents`: Tangents for the edges. Default is `((0,-1),(0,-1))`.
- `arrow_show`: Whether to show arrows on the edges. Default is `false`.
- `arrow_marker`: Marker style for the arrows. Default is `'➤'`.
- `arrow_size`: Size of the arrows. Default is `12`.
- `arrow_shift`: Shift of the arrows. Default is `0.5`.
- `kwargs...`: Additional keyword arguments to pass to the plotting function.

# Returns
- A plot of the network graph with the specified layout and labels.

# Example
```julia
plot_network_coords(eng)
```
"""
function plot_network_coords(
    data::Dict{String,Any};
    show_node_labels=false,
    show_edge_labels=false,
    show_load_labels=false,
    fallback_layout=GraphMakie.Spring(),
    edge_labels_type=:line_id,
    phase="1",
    kwargs...
)

    network_graph, GraphLayout, _ = create_network_graph(data, fallback_layout)

    # Handle labels if required
    nlabels = show_node_labels ? _write_nlabels(network_graph, data) : nothing
    lolabels = show_load_labels ? _write_lolabels(network_graph, data) : nothing
    elabels = show_edge_labels ? edge_labels_type == :line_id ? _write_line_id_elabels(network_graph, data) : _write_results_elabels(network_graph, data, phase) : nothing


    # FORCED NODE FORMATTING:   
    if _is_eng(data)
        _decorate_nodes!(network_graph, data)
        node_color = [props(network_graph, i)[:node_color] for i in 1:nv(network_graph)]
        node_marker = [props(network_graph, i)[:node_marker] for i in 1:nv(network_graph)]
        node_size = [props(network_graph, i)[:marker_size] for i in 1:nv(network_graph)]
        _decorate_edges(network_graph, data)
        edge_color = [get_prop(network_graph, e, :edge_color) for e in edges(network_graph)]
        arrow_show = [get_prop(network_graph, e, :arrow_show) for e in edges(network_graph)]
        arrow_marker = [get_prop(network_graph, e, :arrow_marker) for e in edges(network_graph)]
        arrow_size = [get_prop(network_graph, e, :arrow_size) for e in edges(network_graph)]
        arrow_shift = [get_prop(network_graph, e, :arrow_shift) for e in edges(network_graph)]
    else
        node_color = :black
        node_marker = :circle
        node_size = 10
        edge_color = :black
        arrow_show = false
        arrow_marker = '➤'
        arrow_size = 12
        arrow_shift = 0.5
    end

    # Plot the graph
    return network_graph_plot(
        network_graph;
        layout=GraphLayout, show_node_labels=show_node_labels,
        nlabels=nlabels, #TODO: add extra kwarg for lolabels to re-use nodes
        node_color=node_color,
        node_marker=node_marker,
        node_size=node_size,
        show_edge_labels=show_edge_labels,
        elabels=elabels,
        edge_color=edge_color,
        arrow_show=arrow_show,
        arrow_marker=arrow_marker,
        arrow_size=arrow_size,
        arrow_shift=arrow_shift,
        kwargs...
    )
end

function plot_network_coords!(
    data::Dict{String,Any};
    show_node_labels=false,
    show_edge_labels=false,
    fallback_layout=GraphMakie.Spring(),
    edge_labels_type=:line_id,
    kwargs...
)

    network_graph, GraphLayout, _ = create_network_graph(data, fallback_layout)

    # Handle labels if required
    nlabels = show_node_labels ? _write_nlabels(network_graph, data) : nothing
    elabels = show_edge_labels ? edge_labels_type == :line_id ? _write_line_id_elabels(network_graph, data) : _write_results_elabels(network_graph, data, phase) : nothing

    # FORCED NODE FORMATTING:   
    _decorate_nodes!(network_graph, data)
    node_color = [props(network_graph, i)[:node_color] for i in 1:nv(network_graph)]
    node_marker = [props(network_graph, i)[:node_marker] for i in 1:nv(network_graph)]
    node_size = [props(network_graph, i)[:marker_size] for i in 1:nv(network_graph)]

    _decorate_edges(network_graph, data)
    edge_color = [get_prop(network_graph, e, :edge_color) for e in edges(network_graph)]
    # Plot the graph
    return network_graph_plot!(
        network_graph;
        layout=GraphLayout, show_node_labels=show_node_labels,
        nlabels=nlabels,
        node_color=node_color,
        node_marker=node_marker,
        node_size=node_size,
        show_edge_labels=show_edge_labels,
        elabels=elabels,
        edge_color=edge_color,
        kwargs...
    )
end

"""
    plot_network_map(eng::Dict{String, Any}; show_node_labels=false, show_edge_labels=false, kwargs...)

Plots a network map based on the provided engineering data.

# Arguments
- `eng::Dict{String, Any}`: A dictionary containing the engineering data required to create the network graph.
- `show_node_labels::Bool`: A flag to indicate whether to display labels on the nodes. Default is `false`.
- `show_edge_labels::Bool`: A flag to indicate whether to display labels on the edges. Default is `false`.
- `network_graph::MetaDiGraph`: The network graph to be plotted.
- `GraphLayout::Function`: The layout function for positioning the nodes.
- `tiles_provider`: The tile provider for the map background. Default is `TileProviders.Google(:satelite)`.
- `zoom_lon`: The longitudinal zoom level. Default is `0.0942`.
- `zoom_lat`: The latitudinal zoom level. Default is `0.0942`.
- `makie_backend`: The Makie backend to use for plotting. Default is `WGLMakie`.
- `figure_size`: The size of the figure in pixels. Default is `(1000, 1200)`.
- `nlabels`: Node labels. Default is `nothing`.
- `ilabels`: Internal labels. Default is `nothing`.
- `node_color`: Color of the nodes. Default is `:black`.
- `node_size`: Size of the nodes. Default is `automatic`.
- `node_marker`: Marker style for the nodes. Default is `automatic`.
- `node_strokewidth`: Stroke width of the nodes. Default is `automatic`.
- `show_node_labels`: Whether to show node labels. Default is `false`.
- `elabels`: Edge labels. Default is `nothing`.
- `show_edge_labels`: Whether to show edge labels. Default is `false`.
- `edge_color`: Color of the edges. Default is `:black`.
- `elabels_color`: Color of the edge labels. Default is `:black`.
- `elabels_fontsize`: Font size of the edge labels. Default is `10`.
- `tangents`: Tangents for the edges. Default is `((0,-1),(0,-1))`.
- `arrow_show`: Whether to show arrows on the edges. Default is `false`.
- `arrow_marker`: Marker style for the arrows. Default is `'➤'`.
- `arrow_size`: Size of the arrows. Default is `12`.
- `arrow_shift`: Shift of the arrows. Default is `0.5`.
- `kwargs...`: Additional keyword arguments to customize the plot.

# Returns
- A plot of the network graph with optional node and edge labels.
# Example
```julia
plot_network_map(eng)
`````
"""

function plot_network_map(
    data::Dict{String,Any};
    show_node_labels=false,
    show_edge_labels=false,
    edge_labels_type=:line_id,
    phase="1",
    kwargs...
)
    network_graph, GraphLayout = create_network_graph(data, GraphMakie.Spring())

    # Handle labels if required
    nlabels = show_node_labels ? _write_nlabels(network_graph, data) : nothing
    elabels = show_edge_labels ? edge_labels_type == :line_id ? _write_line_id_elabels(network_graph, data) : _write_results_elabels(network_graph, data, phase) : nothing


    if !isa(GraphLayout, Function)
        @warn "You are attempting to plot a network without coordinates on the map, that is not possible, instead the network tree graph will be plotted"
        return plot_network_coords(data, show_node_labels=show_node_labels, show_edge_labels=show_edge_labels, kwargs...)
    else
        @info "Plotting network map with coordinates on the map -- it is your responsibility to ensure the coordinates are at the correct place"
        if _is_eng(data)
            _decorate_nodes!(network_graph, data)
            node_color = [props(network_graph, i)[:node_color] for i in 1:nv(network_graph)]
            node_marker = [props(network_graph, i)[:node_marker] for i in 1:nv(network_graph)]
            node_size = [props(network_graph, i)[:marker_size] for i in 1:nv(network_graph)]
            _decorate_edges(network_graph, data)
            edge_color = [get_prop(network_graph, e, :edge_color) for e in edges(network_graph)]
            arrow_show = [get_prop(network_graph, e, :arrow_show) for e in edges(network_graph)]
            arrow_marker = [get_prop(network_graph, e, :arrow_marker) for e in edges(network_graph)]
            arrow_size = [get_prop(network_graph, e, :arrow_size) for e in edges(network_graph)]
            arrow_shift = [get_prop(network_graph, e, :arrow_shift) for e in edges(network_graph)]
        else
            node_color = :black
            node_marker = :circle
            node_size = 10
            edge_color = :black
            arrow_show = false
            arrow_marker = '➤'
            arrow_size = 12
            arrow_shift = 0.5
        end
        return network_graph_map_plot(
            network_graph, GraphLayout;
            nlabels=nlabels,
            elabels=elabels,
            show_node_labels=show_node_labels,
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

    end

end

#=
 _______                                                      __      __                     
/       \                                                    /  |    /  |                    
$$$$$$$  |  ______    _______   ______    ______   ______   _$$ |_   $$/   ______   _______  
$$ |  $$ | /      \  /       | /      \  /      \ /      \ / $$   |  /  | /      \ /       \ 
$$ |  $$ |/$$$$$$  |/$$$$$$$/ /$$$$$$  |/$$$$$$  |$$$$$$  |$$$$$$/   $$ |/$$$$$$  |$$$$$$$  |
$$ |  $$ |$$    $$ |$$ |      $$ |  $$ |$$ |  $$/ /    $$ |  $$ | __ $$ |$$ |  $$ |$$ |  $$ |
$$ |__$$ |$$$$$$$$/ $$ \_____ $$ \__$$ |$$ |     /$$$$$$$ |  $$ |/  |$$ |$$ \__$$ |$$ |  $$ |
$$    $$/ $$       |$$       |$$    $$/ $$ |     $$    $$ |  $$  $$/ $$ |$$    $$/ $$ |  $$ |
$$$$$$$/   $$$$$$$/  $$$$$$$/  $$$$$$/  $$/       $$$$$$$/    $$$$/  $$/  $$$$$$/  $$/   $$/ 



=#

function _decorate_nodes!(network_graph::MetaDiGraph, data::Dict{String,Any})
    if _is_eng(data)
        sourcebus = data["voltage_source"]["source"]["bus"]

        for (_, node) in network_graph.vprops
            if node[:bus_id] == Symbol(sourcebus)
                node[:node_color] = :orange
                node[:node_marker] = :star5
                node[:marker_size] = 25
            else
                if !isempty(node[:loads])

                    if length(node[:loads]) == 1
                        if length(node[:loads][1][:connections]) == 1
                            if node[:loads][1][:connections] == [1]
                                node[:node_color] = :red
                                node[:node_marker] = :dtriangle  # `↓`
                            elseif node[:loads][1][:connections] == [2]
                                node[:node_color] = :green
                                node[:node_marker] = :dtriangle  # `↓`
                            elseif node[:loads][1][:connections] == [3]
                                node[:node_color] = :blue
                                node[:node_marker] = :dtriangle  # `↓`
                            elseif node[:loads][1][:connections] == [4]
                                node[:node_color] = :black
                                node[:node_marker] = :dtriangle  # `↓`
                            else
                                error("Unexpected load connections of length $(length(node[:loads])): $(node[:loads][1][:connections]) at node $(string(node[:bus_id]))")
                            end
                        elseif length(node[:loads][1][:connections]) == 2
                            if node[:loads][1][:connections] == [1, 4]
                                node[:node_color] = :red
                                node[:node_marker] = :utriangle  # `↕`
                            elseif node[:loads][1][:connections] == [2, 4]
                                node[:node_color] = :green
                                node[:node_marker] = :utriangle  # `↕`
                            elseif node[:loads][1][:connections] == [3, 4]
                                node[:node_color] = :blue
                                node[:node_marker] = :utriangle  # `↕`
                            elseif node[:loads][1][:connections] == [1, 2]
                                node[:node_color] = :blue
                                node[:node_marker] = :rtriangle  # `↔`
                            elseif node[:loads][1][:connections] == [2, 3]
                                node[:node_color] = :blue
                                node[:node_marker] = :rtriangle  # `↔`
                            elseif node[:loads][1][:connections] == [3, 1]
                                node[:node_color] = :blue
                                node[:node_marker] = :rtriangle  # `↔`
                            else
                                error("Unexpected load connections of length $(length(node[:loads])): $(node[:loads][1][:connections]) at node $(string(node[:bus_id]))")
                            end
                        elseif length(node[:loads][1][:connections]) == 3
                            if node[:loads][1][:connections] == [1, 2, 3]
                                node[:node_color] = :purple
                                node[:node_marker] = :pentagon
                            else
                                error("Unexpected load connections of length $(length(node[:loads])): $(node[:loads][1][:connections]) at node $(string(node[:bus_id]))")
                            end
                        elseif length(node[:loads][1][:connections]) == 4
                            if node[:loads][1][:connections] == [1, 2, 3, 4]
                                node[:node_color] = :purple
                                node[:node_marker] = :pentagon
                            else
                                error("Unexpected load connections of length $(length(node[:loads])): $(node[:loads][1][:connections]) at node $(string(node[:bus_id]))")
                            end
                        else
                            error("Unexpected load connections of length $(length(node[:loads])): $(node[:loads][1][:connections]) at node $(string(node[:bus_id]))")
                        end
                        node[:marker_size] = 10
                    else
                        node[:node_color] = :purple
                        node[:node_marker] = :xcross
                        node[:marker_size] = 10
                    end

                else
                    node[:node_color] = :black
                    node[:node_marker] = :circle
                    node[:marker_size] = 1
                end
            end
        end
    else
        # TODO: MATH related formatting
    end
end


function _decorate_edges(network_graph::MetaDiGraph, data::Dict{String,Any})
    if _is_eng(data)
        for (_, edge) in network_graph.eprops
            # Set default arrow properties
            edge[:arrow_show] = false
            edge[:arrow_marker] = '⠀' #'↡' 
            edge[:arrow_size] = 12
            edge[:arrow_shift] = 0.5

            if !haskey(edge, :t_connections) # I am suing this to determine that it is a transformer not a line
                edge[:edge_color] = :gold
                edge[:arrow_show] = true
                edge[:arrow_marker] = 'Ꝏ'
                edge[:arrow_size] = 20
                edge[:arrow_shift] = 0.5
                
            else

                if !isempty(edge[:t_connections])
                    if length(edge[:t_connections]) == 1
                        if edge[:t_connections] == [1]
                            edge[:edge_color] = :red
                        elseif edge[:t_connections] == [2]
                            edge[:edge_color] = :green
                        elseif edge[:t_connections] == [3]
                            edge[:edge_color] = :blue
                        elseif edge[:t_connections] == [4]
                            edge[:edge_color] = :black
                        else
                            error("Unexpected connections: $(edge[:t_connections])")
                        end
                    elseif length(edge[:t_connections]) == 2
                        if edge[:t_connections] == [1, 4] || edge[:t_connections] == [4, 1]
                            edge[:edge_color] = :red
                        elseif edge[:t_connections] == [2, 4] || edge[:t_connections] == [4, 2]
                            edge[:edge_color] = :green
                        elseif edge[:t_connections] == [3, 4] || edge[:t_connections] == [4, 3]
                            edge[:edge_color] = :blue
                        elseif edge[:t_connections] == [1, 2] || edge[:t_connections] == [2, 1]
                            edge[:edge_color] = :blue
                        elseif edge[:t_connections] == [2, 3] || edge[:t_connections] == [3, 2]
                            edge[:edge_color] = :blue
                        elseif edge[:t_connections] == [3, 1] || edge[:t_connections] == [1, 3]
                            edge[:edge_color] = :blue
                        else
                            error("Unexpected connections: $(edge[:t_connections])")
                        end
                    elseif length(edge[:t_connections]) == 3
                        if edge[:t_connections] == [1, 2, 3]
                            edge[:edge_color] = :purple
                        elseif edge[:t_connections] == [1, 2, 4]
                            edge[:edge_color] = :blue
                        elseif edge[:t_connections] == [2, 3, 4]
                            edge[:edge_color] = :blue
                        elseif edge[:t_connections] == [3, 1, 4]
                            edge[:edge_color] = :blue
                        else
                            error("Unexpected connections: $(edge[:t_connections])")
                        end
                    elseif length(edge[:t_connections]) == 4
                        if length(edge[:t_connections]) == 4 
                            edge[:edge_color] = :purple
                        else
                            error("Unexpected connections: $(edge[:t_connections])")
                        end
                    else
                        edge[:edge_color] = :gray
                    end
                else
                    edge[:edge_color] = :black

                end
            end



        end
    else
        #TODO: MATH related formatting
    end
end
function _decorate_edges(network_graph::MetaDiGraph, data::Dict{String,Any}, phase::String)
    if _is_eng(data)
        for (_, edge) in network_graph.eprops
            if !haskey(edge, :t_connections) # I am suing this to determine that it is a transformer not a line
                edge[:edge_color] = :gold
            else
                if !isempty(edge[:t_connections])
                    if length(edge[:t_connections]) == 1
                        if edge[:t_connections] == [1]
                            edge[:edge_color] = :red
                        elseif edge[:t_connections] == [2]
                            edge[:edge_color] = :white
                        elseif edge[:t_connections] == [3]
                            edge[:edge_color] = :white
                        elseif edge[:t_connections] == [4]
                            edge[:edge_color] = :white
                        else
                            error("Unexpected connections: $(edge[:t_connections])")
                        end
                    elseif length(edge[:t_connections]) == 2
                        if edge[:t_connections] == [1, 4]
                            edge[:edge_color] = :red
                        elseif edge[:t_connections] == [2, 4]
                            edge[:edge_color] = :white
                        elseif edge[:t_connections] == [3, 4]
                            edge[:edge_color] = :white
                        elseif edge[:t_connections] == [1, 2]
                            edge[:edge_color] = :white
                        elseif edge[:t_connections] == [2, 3]
                            edge[:edge_color] = :white
                        elseif edge[:t_connections] == [3, 1]
                            edge[:edge_color] = :white
                        else
                            error("Unexpected connections: $(edge[:t_connections])")
                        end
                    elseif length(edge[:t_connections]) == 3
                        if edge[:t_connections] == [1, 2, 3]
                            edge[:edge_color] = :red
                        elseif edge[:t_connections] == [1, 2, 4]
                            edge[:edge_color] = :red
                        elseif edge[:t_connections] == [2, 3, 4]
                            edge[:edge_color] = :white
                        elseif edge[:t_connections] == [3, 1, 4]
                            edge[:edge_color] = :red
                        else
                            error("Unexpected connections: $(edge[:t_connections])")
                        end
                    elseif length(edge[:t_connections]) == 4
                        if edge[:t_connections] == [1, 2, 3, 4]
                            edge[:edge_color] = :red
                        else
                            error("Unexpected connections: $(edge[:t_connections])")
                        end
                    else
                        edge[:edge_color] = :white
                    end
                else
                    edge[:edge_color] = :white

                end
            end
        end
    else
        #TODO: MATH related formatting
    end
end

#=
__                  __                  __                  __    __                            __  __  __                     
/  |                /  |                /  |                /  |  /  |                          /  |/  |/  |                    
$$ |        ______  $$ |____    ______  $$ |  _______       $$ |  $$ |  ______   _______    ____$$ |$$ |$$/  _______    ______  
$$ |       /      \ $$      \  /      \ $$ | /       |      $$ |__$$ | /      \ /       \  /    $$ |$$ |/  |/       \  /      \ 
$$ |       $$$$$$  |$$$$$$$  |/$$$$$$  |$$ |/$$$$$$$/       $$    $$ | $$$$$$  |$$$$$$$  |/$$$$$$$ |$$ |$$ |$$$$$$$  |/$$$$$$  |
$$ |       /    $$ |$$ |  $$ |$$    $$ |$$ |$$      \       $$$$$$$$ | /    $$ |$$ |  $$ |$$ |  $$ |$$ |$$ |$$ |  $$ |$$ |  $$ |
$$ |_____ /$$$$$$$ |$$ |__$$ |$$$$$$$$/ $$ | $$$$$$  |      $$ |  $$ |/$$$$$$$ |$$ |  $$ |$$ \__$$ |$$ |$$ |$$ |  $$ |$$ \__$$ |
$$       |$$    $$ |$$    $$/ $$       |$$ |/     $$/       $$ |  $$ |$$    $$ |$$ |  $$ |$$    $$ |$$ |$$ |$$ |  $$ |$$    $$ |
$$$$$$$$/  $$$$$$$/ $$$$$$$/   $$$$$$$/ $$/ $$$$$$$/        $$/   $$/  $$$$$$$/ $$/   $$/  $$$$$$$/ $$/ $$/ $$/   $$/  $$$$$$$ |
                                                                                                                      /  \__$$ |
                                                                                                                      $$    $$/                                                                                                                      $$$$$$/  
=#

## NODES
function _write_nlabels(network_graph, data)
    nlabels = []
    for i in 1:nv(network_graph)
        push!(nlabels,
            string(network_graph[i, :bus_id])
        )
    end
    return nlabels
end

## LOADS

function _write_lolabels(network_graph, data)
    lolabels = []
    for i in 1:nv(network_graph)
        if !isempty(network_graph.vprops[i][:loads])
            load_str = ""
            for load in network_graph.vprops[i][:loads]
                load_str *= "$(load[:bus])\n"
            end
            push!(lolabels, load_str)
        else
            push!(lolabels, "")
        end
    end
    @show lolabels
    return lolabels
end


## EDGES


function _write_line_id_elabels(network_graph, data)
    if _is_eng(data)
        [string(get_prop(network_graph, e, :line_id)) for e in edges(network_graph)]
    else
        [string(get_prop(network_graph, e, :branch_id)) for e in edges(network_graph)]
    end
end


function _check_bases(eng_res::Dict{String,Any})
    if haskey(eng_res, "bases")
        return eng_res["bases"]["is_perunit"], eng_res["bases"]["vbase_V"], eng_res["bases"]["sbase_VA"], eng_res["bases"]["Zbase_Ω"], eng_res["bases"]["Ibase_A"]
    else
        return false, 0.0, 0.0, 0.0, 0.0
    end

end


function _write_results_elabels(network_graph, data, phase)
    if !_is_eng(data)
        error("The data is not an engineering model, also it has to have the results to be able to plot them the edge labels")
    end

    try
        elabels = []
        _, vbase_V, _, _, Ibase_A = _check_bases(data)

        for e in edges(network_graph)

            var_i_dict = get_prop(network_graph, e, :I_f) # variable current
            var_v_dict = get_prop(network_graph, e, :V_f) # variable voltage

            if haskey(var_i_dict, Symbol(phase)) && haskey(var_v_dict, Symbol(phase))
                push!(
                    elabels,
                    "V$phase =$(round(abs.(var_v_dict[Symbol(phase)])*vbase_V, digits=2)) ∠ $(round(rad2deg.(angle.(var_v_dict[Symbol(phase)])), digits=2))\n
     I$phase=$(round(abs.(var_i_dict[Symbol(phase)])*Ibase_A, digits=2)) ∠ $(round(rad2deg.(angle.(var_i_dict[Symbol(phase)])), digits=2))"
                )
            else
                push!(elabels, "N/A")
            end
        end

        return elabels

    catch e
        error("maybe the data model you passed doesn't have the results?  You need to run for example `eng_res, math, PF = pf_results(eng);` and then pass the eng_res which has results \n error is: \n $e")
    end
end









# Re-export plotting functions from parent module
export create_graph
export plot_network_tree, plot_network_tree!
export plot_network_coords, plot_network_coords!
export plot_network_map
export bus_phasor, bus_phasor!

end # module PMDGraph
