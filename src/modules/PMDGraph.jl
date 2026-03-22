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


function _is_eng_graph(graph::AbstractMetaGraph)

    if haskey(first(graph.eprops).second, :branch_id)
        return false
    else
        return true
    end
    error("I don't know if this graph is for a MATHEMATICAL or ENGINEERING model. I usually check if it has `:branch_id` or `:line_id` in the edge properties to tell which one it is.")
end


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
"""
function create_network_graph_eng(eng::Dict{String,Any}, fallback_layout)
    eng_sym = convert_keys_to_symbols(deepcopy(eng))
    network_graph = MetaDiGraph()

    # --- 1. Add Vertices (Buses) ---
    # Determine Root Bus (Source Bus)
    ref_bus_name = Symbol("")
    if haskey(eng_sym, :voltage_source) && haskey(eng_sym[:voltage_source], :source)
        ref_bus_name = Symbol(eng_sym[:voltage_source][:source][:bus])
    end

    # Add root first if it exists
    if ref_bus_name != Symbol("") && haskey(eng_sym[:bus], ref_bus_name)
        bus = eng_sym[:bus][ref_bus_name]
        bus[:bus_id] = ref_bus_name
        bus[:loads] = []
        bus[:shunts] = []
        bus[:voltage_sources] = []
        add_vertex!(network_graph, bus)
    end

    for (b, bus) in eng_sym[:bus]
        # Skip if already added
        if b == ref_bus_name
            continue
        end

        # Ensure bus_id is stored consistently
        bus[:bus_id] = b

        # Initialize connection lists
        bus[:loads] = []
        bus[:shunts] = []
        bus[:voltage_sources] = []

        add_vertex!(network_graph, bus)
    end

    # Use bus_id as the indexing property immediately after adding all vertices
    set_indexing_prop!(network_graph, :bus_id)

    # --- 2. Enrich Bus Data ---
    # Attach loads
    if haskey(eng_sym, :load)
        for (load_id, load) in eng_sym[:load]
            load[:load_id] = load_id
            bus_id = Symbol(load[:bus])
            try
                vertex_idx = network_graph[bus_id, :bus_id]
                push!(network_graph.vprops[vertex_idx][:loads], load)
            catch
                @warn "Load connected to missing bus: $bus_id"
            end
        end
    end

    # Attach shunts
    if haskey(eng_sym, :shunt)
        for (_, shunt) in eng_sym[:shunt]
            bus_id = Symbol(shunt[:bus])
            try
                vertex_idx = network_graph[bus_id, :bus_id]
                push!(network_graph.vprops[vertex_idx][:shunts], shunt)
            catch
                @warn "Shunt connected to missing bus: $bus_id"
            end
        end
    end

    # Attach voltage sources
    if haskey(eng_sym, :voltage_source)
        for (_, vs) in eng_sym[:voltage_source]
            bus_id = Symbol(vs[:bus])
            try
                vertex_idx = network_graph[bus_id, :bus_id]
                push!(network_graph.vprops[vertex_idx][:voltage_sources], vs)
            catch
                @warn "Voltage source connected to missing bus: $bus_id"
            end
        end
    end

    # --- 3. Add Edges (Lines & Transformers) ---
    # Used to track connectivity for coordinate inference
    adj_list = Dict{Symbol,Vector{Symbol}}()

    # Helper to add edge and track adjacency
    function safe_add_edge!(f_id, t_id, data)
        try
            f_v = network_graph[f_id, :bus_id]
            t_v = network_graph[t_id, :bus_id]
            add_edge!(network_graph, f_v, t_v, data)

            # Track adjacency for coordinate back-propagation
            push!(get!(adj_list, f_id, []), t_id)
            push!(get!(adj_list, t_id, []), f_id)
        catch e
            @warn "Could not add edge between $f_id and $t_id: $e"
        end
    end

    for (l, line) in eng_sym[:line]
        line[:line_id] = l
        if haskey(line, :linecode) && haskey(eng_sym, :linecode)
            line[:linecodes] = eng_sym[:linecode][Symbol(line[:linecode])]
        end
        safe_add_edge!(Symbol(line[:f_bus]), Symbol(line[:t_bus]), line)
    end

    if haskey(eng_sym, :transformer)
        for (t_id, trans) in eng_sym[:transformer]
            trans[:line_id] = t_id
            # Transformers in eng model connect array of buses
            buses = trans[:bus]
            if length(buses) >= 2
                safe_add_edge!(Symbol(buses[1]), Symbol(buses[2]), trans)
            end
        end
    end

    if haskey(eng_sym, :switch)
        for (sw_id, sw) in eng_sym[:switch]
            sw[:line_id] = sw_id
            sw[:is_switch] = true
            safe_add_edge!(Symbol(sw[:f_bus]), Symbol(sw[:t_bus]), sw)
        end
    end

    # --- 4. Coordinate Handling & Root Fix ---

    # Collect existing coordinates
    # We iterate over vertices to respect the graph order (1..nv) for layout vectors
    layouting_vector = Vector{Tuple{Float64,Float64}}(undef, nv(network_graph))
    has_coords = fill(false, nv(network_graph))

    for v in vertices(network_graph)
        bus_props = props(network_graph, v)
        if haskey(bus_props, :lon) && haskey(bus_props, :lat)
            layouting_vector[v] = (bus_props[:lon], bus_props[:lat])
            has_coords[v] = true
        end
    end

    # Fix missing coordinates for root buses (or any bus) by looking at neighbors
    # This is a simple 1-hop look-ahead. Can be expanded to BFS if needed.

    # Identify likely root buses (those with voltage sources)
    root_candidates = []
    if haskey(eng_sym, :voltage_source)
        for (_, vs) in eng_sym[:voltage_source]
            push!(root_candidates, Symbol(vs[:bus]))
        end
    end

    # Also iterate all buses to patch holes generally
    for _ in 1:2 # Two passes for propagation
        for v in vertices(network_graph)
            if !has_coords[v]
                # Look for a neighbor with coordinates
                v_sym = get_prop(network_graph, v, :bus_id)
                neighbors = get(adj_list, v_sym, [])
                for neighbor_id in neighbors
                    try
                        v_neighbor = network_graph[neighbor_id, :bus_id]
                        if has_coords[v_neighbor]
                            layouting_vector[v] = layouting_vector[v_neighbor]
                            has_coords[v] = true
                            # Update the property in the graph as well for consistency
                            set_prop!(network_graph, v, :lon, layouting_vector[v_neighbor][1])
                            set_prop!(network_graph, v, :lat, layouting_vector[v_neighbor][2])
                            break
                        end
                    catch
                    end
                end
            end
        end
    end

    # Layout function construction
    if sum(has_coords) > 0 # At least some coordinates exist
        GraphLayout = function (g)
            final_layout = []
            for i in 1:nv(g)
                if has_coords[i]
                    push!(final_layout, layouting_vector[i])
                else
                    push!(final_layout, (0.0, 0.0)) # Placeholder
                end
            end
            return final_layout
        end
        if sum(has_coords) < nv(network_graph)
            @warn "Only $(sum(has_coords)) of $(nv(network_graph)) buses have coordinates. Use fallback for missing ones?"
        end
    else
        @warn "No coordinates found in ENGINEERING model. Using fallback layout."
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
"""
function create_network_graph_math(math::Dict{String,Any}, fallback_layout)
    math_sym = convert_keys_to_symbols(deepcopy(math))
    network_graph = MetaDiGraph()

    # --- 1. Add Vertices (Buses) ---
    # Prioritize adding Root/Virtual buses first for layout algorithms (Buchheim)

    # 1. Find physical reference bus (type 3)
    phys_ref_id = nothing
    for (b, bus) in math_sym[:bus]
        if get(bus, :bus_type, 1) == 3
            phys_ref_id = b
            break
        end
    end

    # 2. Check for "virtual bus" upstream (connected via branch where ref is t_bus)
    virt_bus_id = nothing
    if phys_ref_id !== nothing
        for (_, branch) in math_sym[:branch]
            # math keys might differ in type (int vs symbol), check carefully or rely on consistency
            if branch[:t_bus] == phys_ref_id
                virt_bus_id = branch[:f_bus]
                break
            end
        end
    end

    priority_buses = Set()

    # helper to add bus
    function add_bus_node!(b_id)
        if haskey(math_sym[:bus], b_id) && !(b_id in priority_buses)
            bus = math_sym[:bus][b_id]
            bus[:bus_id] = b_id
            bus[:loads] = []
            bus[:shunts] = []
            bus[:gens] = []
            add_vertex!(network_graph, bus)
            push!(priority_buses, b_id)
        end
    end

    # Add in order: Virtual -> Physical Ref
    if virt_bus_id !== nothing
        add_bus_node!(virt_bus_id)
    end
    if phys_ref_id !== nothing
        add_bus_node!(phys_ref_id)
    end

    for (b, bus) in math_sym[:bus]
        if !(b in priority_buses)
            bus[:bus_id] = b
            bus[:loads] = []
            bus[:shunts] = []
            bus[:gens] = []

            add_vertex!(network_graph, bus)
        end
    end

    # Set indexing property
    set_indexing_prop!(network_graph, :bus_id)

    # --- 2. Enrich Bus Data ---
    if haskey(math_sym, :load)
        for (load_id, load) in math_sym[:load]
            try
                load[:load_id] = load_id
                bus_id = Symbol(load[:load_bus]) # math model often uses ints, ensure key match

                # Robust lookup: try direct, then Symbol, then Int
                v_idx = nothing
                try
                    v_idx = network_graph[bus_id, :bus_id]
                catch
                end

                if !isnothing(v_idx)
                    push!(network_graph.vprops[v_idx][:loads], load)
                end
            catch e
                @debug "Failed to attach load: $e"
            end
        end
    end

    if haskey(math_sym, :shunt)
        for (_, shunt) in math_sym[:shunt]
            try
                bus_id = Symbol(shunt[:shunt_bus])
                v_idx = nothing
                try
                    v_idx = network_graph[bus_id, :bus_id]
                catch
                end
                if !isnothing(v_idx)
                    push!(network_graph.vprops[v_idx][:shunts], shunt)
                end
            catch
            end
        end
    end

    if haskey(math_sym, :gen)
        for (_, gen) in math_sym[:gen]
            try
                bus_id = Symbol(gen[:gen_bus])
                v_idx = nothing
                try
                    v_idx = network_graph[bus_id, :bus_id]
                catch
                end
                if !isnothing(v_idx)
                    push!(network_graph.vprops[v_idx][:gens], gen)
                end
            catch
            end
        end
    end

    # --- 3. Add Edges (Branches) ---
    adj_list = Dict{Any,Vector{Any}}()

    for (l, branch) in math_sym[:branch]
        branch[:branch_id] = l

        f_bus = Symbol(branch[:f_bus])
        t_bus = Symbol(branch[:t_bus])

        if startswith(string(get(branch, :name, "")), "_virtual_branch") && endswith(string(get(branch, :name, "")), "_2")
            @debug "Branch $(l) appears to be a reversed transformer branch (ends with '_2'), flipping direction for graph edge"
            f_bus, t_bus = t_bus, f_bus
        end

        try
            f_vertex = network_graph[f_bus, :bus_id]
            t_vertex = network_graph[t_bus, :bus_id]

            add_edge!(network_graph, f_vertex, t_vertex, branch)

            push!(get!(adj_list, f_bus, []), t_bus)
            push!(get!(adj_list, t_bus, []), f_bus)
        catch e
            @warn "Error adding edge for branch $(l) ($f_bus -> $t_bus): $e"
        end
    end

    if haskey(math_sym, :transformer)
        for (t, trans) in math_sym[:transformer]
            trans[:branch_id] = t
            trans[:is_transformer] = true

            f_bus = Symbol(trans[:f_bus])
            t_bus = Symbol(trans[:t_bus])

            if startswith(string(get(trans, :name, "")), "_virtual_transformer") && endswith(string(get(trans, :name, "")), ".2")
                @debug "Transformer $(t) appears to be a reversed branch (ends with '.2'), flipping direction for graph edge"
                f_bus, t_bus = t_bus, f_bus
            end

            try
                f_vertex = network_graph[f_bus, :bus_id]
                t_vertex = network_graph[t_bus, :bus_id]

                if has_edge(network_graph, f_vertex, t_vertex)
                    for (k, v) in trans
                        set_prop!(network_graph, f_vertex, t_vertex, k, v)
                    end
                else
                    add_edge!(network_graph, f_vertex, t_vertex, trans)
                    push!(get!(adj_list, f_bus, []), t_bus)
                    push!(get!(adj_list, t_bus, []), f_bus)
                end
            catch e
                @warn "Error adding edge for transformer $(t) ($f_bus -> $t_bus): $e"
            end
        end
    end

    if haskey(math_sym, :switch)
        for (sw_id, sw) in math_sym[:switch]
            sw[:branch_id] = sw_id
            sw[:is_switch] = true

            f_bus = Symbol(sw[:f_bus])
            t_bus = Symbol(sw[:t_bus])

            try
                f_vertex = network_graph[f_bus, :bus_id]
                t_vertex = network_graph[t_bus, :bus_id]

                if has_edge(network_graph, f_vertex, t_vertex)
                    for (k, v) in sw
                        set_prop!(network_graph, f_vertex, t_vertex, k, v)
                    end
                else
                    add_edge!(network_graph, f_vertex, t_vertex, sw)
                    push!(get!(adj_list, f_bus, []), t_bus)
                    push!(get!(adj_list, t_bus, []), f_bus)
                end
            catch e
                @warn "Error adding edge for switch $(sw_id) ($f_bus -> $t_bus): $e"
            end
        end
    end

    # --- 4. Coordinate Handling ---
    layouting_vector = Vector{Tuple{Float64,Float64}}(undef, nv(network_graph))
    has_coords = fill(false, nv(network_graph))

    for v in vertices(network_graph)
        bus_props = props(network_graph, v)
        if haskey(bus_props, :lon) && haskey(bus_props, :lat)
            layouting_vector[v] = (bus_props[:lon], bus_props[:lat])
            has_coords[v] = true
        end
    end

    # Propagate coordinates to any node missing them from neighbors
    # (Covers RefBus -> VirtualGenBus and similar cases)
    for _ in 1:2
        for v in vertices(network_graph)
            if !has_coords[v]
                # Check neighbors
                bus_id = get_prop(network_graph, v, :bus_id)
                neighbors = get(adj_list, bus_id, [])

                for n_id in neighbors
                    try
                        n_v = network_graph[n_id, :bus_id]
                        if has_coords[n_v]
                            layouting_vector[v] = layouting_vector[n_v]
                            has_coords[v] = true
                            set_prop!(network_graph, v, :lon, layouting_vector[n_v][1])
                            set_prop!(network_graph, v, :lat, layouting_vector[n_v][2])
                            break
                        end
                    catch
                    end
                end
            end
        end
    end

    if sum(has_coords) > 0
        GraphLayout = function (g)
            final_layout = []
            for i in 1:nv(g)
                if has_coords[i]
                    push!(final_layout, layouting_vector[i])
                else
                    push!(final_layout, (0.0, 0.0))
                end
            end
            return final_layout
        end
    else
        @warn "No coordinates found in MATHEMATICAL model. Using fallback layout."
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
    if _is_eng_graph(G)
        Eid = findall(e -> e[:line_id] == Symbol(edge_id), G.eprops)
    else
        Eid = findall(e -> e[:branch_id] == Symbol(edge_id), G.eprops)
    end
    return G.eprops[Eid...]
end
function get_graph_edge(G, edge_id, key)
    if _is_eng_graph(G)
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
    smart_layout(g)

A smart layout function that checks if the graph `g` allows for a tree-based layout.
If `g` is a tree (weakly connected and |E| = |V| - 1), it attempts to use `GraphMakie.Buchheim()`.
If `Buchheim` fails or if the graph is not a tree, it falls back to `GraphMakie.Stress(Ptype=Float32)()`.
"""
function smart_layout(g)
    is_tree = false
    try
        if nv(g) > 0 && ne(g) == nv(g) - 1
            if is_weakly_connected(g)
                is_tree = true
            end
        end
    catch e
        @warn "Tree check failed: $e"
    end

    if is_tree
        try
            @debug "Graph appears to be a tree, using Buchheim layout"
            return GraphMakie.Buchheim()(g)
        catch e
            @warn "Buchheim layout failed, falling back to Stress" exception = e
            return GraphMakie.Stress(Ptype=Float32)(g)
        end
    else
        return GraphMakie.Stress(Ptype=Float32)(g)
    end
end




"""
    network_graph_plot(network_graph::MetaDiGraph; layout=smart_layout, figure_size=(1000, 1200), node_color=:black, node_size=automatic, node_marker=automatic, node_strokewidth=automatic, show_node_labels=false, show_edge_labels=false, edge_color=:black, elabels_color=:black, elabels_fontsize=10, tangents=((0,-1),(0,-1)), arrow_show=false, arrow_marker='➤', arrow_size=50, arrow_shift=0.5)

Plots the given network graph using the specified layout and styling options.

# Arguments
- `network_graph::MetaDiGraph`: The network graph to be plotted.

# Keyword Arguments
- `layout`: The layout algorithm to use for positioning the nodes (default: `smart_layout`, which uses `Buchheim` for trees and `Spring` otherwise).
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
    layout=smart_layout, #layout=GraphMakie.Spring(),
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
    layout=smart_layout,
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
    layout=smart_layout,
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
    show_load_labels=false,
    layout=smart_layout,
    edge_labels_type=:line_id,
    phase="1",
    kwargs...
)
    # Create the network meta graph 
    network_graph, _, _ = create_network_graph(data, layout)

    # Handle labels if required
    nlabels = show_node_labels ? _write_nlabels(network_graph, data) : nothing
    lolabels = show_load_labels ? _write_lolabels(network_graph, data) : nothing
    elabels = show_edge_labels ? edge_labels_type == :line_id ? _write_line_id_elabels(network_graph, data) : _write_results_elabels(network_graph, data, phase) : nothing

    if show_load_labels
        if show_node_labels
            nlabels = [string(n, "\n", l) for (n, l) in zip(nlabels, lolabels)]
        else
            nlabels = lolabels
        end
    end

    # FORCED NODE FORMATTING:
    _decorate_nodes!(network_graph, data)
    node_color = [props(network_graph, i)[:node_color] for i in 1:nv(network_graph)]
    node_marker = [props(network_graph, i)[:node_marker] for i in 1:nv(network_graph)]
    node_size = [props(network_graph, i)[:marker_size] for i in 1:nv(network_graph)]

    _decorate_edges!(network_graph, data)
    edge_color = [get_prop(network_graph, e, :edge_color) for e in edges(network_graph)]
    # Calculating arrow_show as true (globally enabled) but controlling visibility via arrow_size=0
    # This avoids "non-boolean (Vector{Bool}) used in boolean context" error in Makie
    arrow_show = true
    # arrow_show = [get_prop(network_graph, e, :arrow_show) for e in edges(network_graph)]

    arrow_marker = [get_prop(network_graph, e, :arrow_marker) for e in edges(network_graph)]
    arrow_size = [get_prop(network_graph, e, :arrow_size) for e in edges(network_graph)]
    arrow_shift = [get_prop(network_graph, e, :arrow_shift) for e in edges(network_graph)]
    # plot and return the network 
    return network_graph_plot(
        network_graph;
        layout=layout,
        figure_size=figure_size,
        show_node_labels=show_node_labels || show_load_labels,
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
    show_load_labels=false,
    layout=smart_layout,
    edge_labels_type=:line_id,
    phase="1",
    kwargs...
)
    # Create the network meta graph 
    network_graph, _, _ = create_network_graph(data, layout)

    # Handle labels if required
    nlabels = show_node_labels ? _write_nlabels(network_graph, data) : nothing
    lolabels = show_load_labels ? _write_lolabels(network_graph, data) : nothing
    elabels = show_edge_labels ? edge_labels_type == :line_id ? _write_line_id_elabels(network_graph, data) : _write_results_elabels(network_graph, data, phase) : nothing

    if show_load_labels
        if show_node_labels
            nlabels = [string(n, "\n", l) for (n, l) in zip(nlabels, lolabels)]
        else
            nlabels = lolabels
        end
    end

    # FORCED NODE FORMATTING:
    _decorate_nodes!(network_graph, data)
    node_color = [props(network_graph, i)[:node_color] for i in 1:nv(network_graph)]
    node_marker = [props(network_graph, i)[:node_marker] for i in 1:nv(network_graph)]
    node_size = [props(network_graph, i)[:marker_size] for i in 1:nv(network_graph)]

    _decorate_edges!(network_graph, data)
    edge_color = [get_prop(network_graph, e, :edge_color) for e in edges(network_graph)]
    # Calculating arrow_show as true (globally enabled) but controlling visibility via arrow_size=0
    # This avoids "non-boolean (Vector{Bool}) used in boolean context" error in Makie
    arrow_show = true
    # arrow_show = [get_prop(network_graph, e, :arrow_show) for e in edges(network_graph)]

    arrow_marker = [get_prop(network_graph, e, :arrow_marker) for e in edges(network_graph)]
    arrow_size = [get_prop(network_graph, e, :arrow_size) for e in edges(network_graph)]
    arrow_shift = [get_prop(network_graph, e, :arrow_shift) for e in edges(network_graph)]
    # plot and return the network 
    return network_graph_plot!(
        ax,
        network_graph;
        layout=layout,
        figure_size=figure_size,
        show_node_labels=show_node_labels || show_load_labels,
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
    plot_network_coords(eng::Dict{String,Any}; show_node_labels=false, show_edge_labels=false, fallback_layout=smart_layout, kwargs...)

Plot a network graph with optional node and edge labels.

# Arguments
- `eng::Dict{String,Any}`: The engineering data used to create the network graph.
- `show_node_labels::Bool`: Whether to display labels for the nodes. Default is `false`.
- `show_edge_labels::Bool`: Whether to display labels for the edges. Default is `false`.
- `fallback_layout`: The layout algorithm to use if no specific layout is provided. Default is `smart_layout`.
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
    fallback_layout=smart_layout,
    edge_labels_type=:line_id,
    phase="1",
    kwargs...
)

    network_graph, GraphLayout, _ = create_network_graph(data, fallback_layout)

    # Handle labels if required
    nlabels = show_node_labels ? _write_nlabels(network_graph, data) : nothing
    lolabels = show_load_labels ? _write_lolabels(network_graph, data) : nothing
    elabels = show_edge_labels ? edge_labels_type == :line_id ? _write_line_id_elabels(network_graph, data) : _write_results_elabels(network_graph, data, phase) : nothing

    if show_load_labels
        if show_node_labels
            nlabels = [string(n, "\n", l) for (n, l) in zip(nlabels, lolabels)]
        else
            nlabels = lolabels
        end
    end


    # FORCED NODE FORMATTING:   
    _decorate_nodes!(network_graph, data)
    node_color = [props(network_graph, i)[:node_color] for i in 1:nv(network_graph)]
    node_marker = [props(network_graph, i)[:node_marker] for i in 1:nv(network_graph)]
    node_size = [props(network_graph, i)[:marker_size] for i in 1:nv(network_graph)]

    _decorate_edges!(network_graph, data)
    edge_color = [get_prop(network_graph, e, :edge_color) for e in edges(network_graph)]
    # Calculating arrow_show as true (globally enabled) but controlling visibility via arrow_size=0
    # This avoids "non-boolean (Vector{Bool}) used in boolean context" error in Makie
    arrow_show = true
    # arrow_show = [get_prop(network_graph, e, :arrow_show) for e in edges(network_graph)]

    arrow_marker = [get_prop(network_graph, e, :arrow_marker) for e in edges(network_graph)]
    arrow_size = [get_prop(network_graph, e, :arrow_size) for e in edges(network_graph)]
    arrow_shift = [get_prop(network_graph, e, :arrow_shift) for e in edges(network_graph)]

    # Plot the graph
    return network_graph_plot(
        network_graph;
        layout=GraphLayout, show_node_labels=show_node_labels || show_load_labels,
        nlabels=nlabels,
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
    fallback_layout=smart_layout,
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

    _decorate_edges!(network_graph, data)
    edge_color = [get_prop(network_graph, e, :edge_color) for e in edges(network_graph)]
    arrow_show = [get_prop(network_graph, e, :arrow_show) for e in edges(network_graph)]
    arrow_marker = [get_prop(network_graph, e, :arrow_marker) for e in edges(network_graph)]
    arrow_size = [get_prop(network_graph, e, :arrow_size) for e in edges(network_graph)]
    arrow_shift = [get_prop(network_graph, e, :arrow_shift) for e in edges(network_graph)]

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
        arrow_show=arrow_show,
        arrow_marker=arrow_marker,
        arrow_size=arrow_size,
        arrow_shift=arrow_shift,
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
```
"""
function plot_network_map(
    data::Dict{String,Any};
    show_node_labels=false,
    show_edge_labels=false,
    show_load_labels=false,
    edge_labels_type=:line_id,
    phase="1",
    kwargs...
)
    network_graph, GraphLayout = create_network_graph(data, smart_layout)

    # Handle labels if required
    nlabels = show_node_labels ? _write_nlabels(network_graph, data) : nothing
    lolabels = show_load_labels ? _write_lolabels(network_graph, data) : nothing
    elabels = show_edge_labels ? edge_labels_type == :line_id ? _write_line_id_elabels(network_graph, data) : _write_results_elabels(network_graph, data, phase) : nothing

    if show_load_labels
        if show_node_labels
            nlabels = [string(n, "\n", l) for (n, l) in zip(nlabels, lolabels)]
        else
            nlabels = lolabels
        end
    end


    if !isa(GraphLayout, Function)
        @warn "You are attempting to plot a network without coordinates on the map, that is not possible, instead the network tree graph will be plotted"
        return plot_network_coords(data, show_node_labels=show_node_labels, show_edge_labels=show_edge_labels, kwargs...)
    else
        @info "Plotting network map with coordinates on the map -- it is your responsibility to ensure the coordinates are at the correct place"
        _decorate_nodes!(network_graph, data)
        node_color = [props(network_graph, i)[:node_color] for i in 1:nv(network_graph)]
        node_marker = [props(network_graph, i)[:node_marker] for i in 1:nv(network_graph)]
        node_size = [props(network_graph, i)[:marker_size] for i in 1:nv(network_graph)]

        _decorate_edges!(network_graph, data)
        edge_color = [get_prop(network_graph, e, :edge_color) for e in edges(network_graph)]
        # Calculating arrow_show as true (globally enabled) but controlling visibility via arrow_size=0
        # This avoids "non-boolean (Vector{Bool}) used in boolean context" error in Makie
        arrow_show = true
        # arrow_show = [get_prop(network_graph, e, :arrow_show) for e in edges(network_graph)]

        arrow_marker = [get_prop(network_graph, e, :arrow_marker) for e in edges(network_graph)]
        arrow_size = [get_prop(network_graph, e, :arrow_size) for e in edges(network_graph)]
        arrow_shift = [get_prop(network_graph, e, :arrow_shift) for e in edges(network_graph)]
        return network_graph_map_plot(
            network_graph, GraphLayout;
            nlabels=nlabels,
            elabels=elabels,
            show_node_labels=show_node_labels || show_load_labels,
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
                                Pliers.warning_text("Unexpected load connections of length $(length(node[:loads])): $(node[:loads][1][:connections]) at node $(string(node[:bus_id]))")
                                node[:node_color] = :gray
                                node[:node_marker] = '?'
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
                                Pliers.warning_text("Unexpected load connections of length $(length(node[:loads])): $(node[:loads][1][:connections]) at node $(string(node[:bus_id]))")
                                node[:node_color] = :gray
                                node[:node_marker] = '?'
                            end
                        elseif length(node[:loads][1][:connections]) == 3
                            if node[:loads][1][:connections] == [1, 2, 3]
                                node[:node_color] = :purple
                                node[:node_marker] = :pentagon
                            else
                                Pliers.warning_text("Unexpected load connections of length $(length(node[:loads])): $(node[:loads][1][:connections]) at node $(string(node[:bus_id]))")
                                node[:node_color] = :gray
                                node[:node_marker] = '?'
                            end
                        elseif length(node[:loads][1][:connections]) == 4
                            if node[:loads][1][:connections] == [1, 2, 3, 4]
                                node[:node_color] = :purple
                                node[:node_marker] = :pentagon
                            else
                                Pliers.warning_text("Unexpected load connections of length $(length(node[:loads])): $(node[:loads][1][:connections]) at node $(string(node[:bus_id]))")
                                node[:node_color] = :gray
                                node[:node_marker] = '?'
                            end
                        else
                            Pliers.warning_text("Unexpected load connections of length $(length(node[:loads])): $(node[:loads][1][:connections]) at node $(string(node[:bus_id]))")
                            node[:node_color] = :gray
                            node[:node_marker] = '?'
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
        # MATH related formatting
        for (_, node) in network_graph.vprops
            # Reference Bus (Type 3)
            # Check if bus_type is 3 (Swing/Slack)
            is_ref = false
            if haskey(node, :bus_type) && node[:bus_type] == 3
                is_ref = true
            end

            if is_ref
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
                                Pliers.warning_text("Unexpected load connections of length $(length(node[:loads])): $(node[:loads][1][:connections]) at node $(string(node[:bus_id]))")
                                node[:node_color] = :gray
                                node[:node_marker] = '?'
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
                                Pliers.warning_text("Unexpected load connections of length $(length(node[:loads])): $(node[:loads][1][:connections]) at node $(string(node[:bus_id]))")
                                node[:node_color] = :gray
                                node[:node_marker] = '?'
                            end
                        elseif length(node[:loads][1][:connections]) == 3
                            if node[:loads][1][:connections] == [1, 2, 3]
                                node[:node_color] = :purple
                                node[:node_marker] = :pentagon
                            else
                                Pliers.warning_text("Unexpected load connections of length $(length(node[:loads])): $(node[:loads][1][:connections]) at node $(string(node[:bus_id]))")
                                node[:node_color] = :gray
                                node[:node_marker] = '?'
                            end
                        elseif length(node[:loads][1][:connections]) == 4
                            if node[:loads][1][:connections] == [1, 2, 3, 4]
                                node[:node_color] = :purple
                                node[:node_marker] = :pentagon
                            else
                                Pliers.warning_text("Unexpected load connections of length $(length(node[:loads])): $(node[:loads][1][:connections]) at node $(string(node[:bus_id]))")
                                node[:node_color] = :gray
                                node[:node_marker] = '?'
                            end
                        else
                            Pliers.warning_text("Unexpected load connections of length $(length(node[:loads])): $(node[:loads][1][:connections]) at node $(string(node[:bus_id]))")
                            node[:node_color] = :gray
                            node[:node_marker] = '?'
                        end
                        node[:marker_size] = 10
                    else
                        node[:node_color] = :purple
                        node[:node_marker] = :xcross
                        node[:marker_size] = 10
                    end
                else
                    node[:node_color] = :gray
                    node[:node_marker] = :circle
                    node[:marker_size] = 1
                end
            end
        end
    end
end


function _decorate_edges!(network_graph::MetaDiGraph, data::Dict{String,Any})
    if _is_eng(data)
        for (_, edge) in network_graph.eprops
            # Set default arrow properties
            edge[:arrow_show] = false
            edge[:arrow_marker] = '⠀'  # '→' #'↡' 
            edge[:arrow_size] = 0
            edge[:arrow_shift] = 0.5

            if haskey(edge, :is_switch) && edge[:is_switch]
                # Switch decoration
                edge[:edge_color] = :teal
                edge[:arrow_show] = true
                edge[:arrow_marker] = string(edge[:state]) == "OPEN" ? '□' : '■'
                edge[:arrow_size] = 20
                edge[:arrow_shift] = 0.5

            elseif !haskey(edge, :t_connections) # I am using this to determine that it is a transformer not a line
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
                            Pliers.warning_text("Unexpected connections: $(edge[:t_connections])")
                            edge[:edge_color] = :gray
                        end
                    elseif length(edge[:t_connections]) == 2
                        if edge[:t_connections] == [1, 4] || edge[:t_connections] == [4, 1]
                            edge[:edge_color] = :red
                        elseif edge[:t_connections] == [2, 4] || edge[:t_connections] == [4, 2]
                            edge[:edge_color] = :green
                        elseif edge[:t_connections] == [3, 4] || edge[:t_connections] == [4, 3]
                            edge[:edge_color] = :blue
                        elseif edge[:t_connections] == [1, 2] || edge[:t_connections] == [2, 1]
                            edge[:edge_color] = :cyan
                        elseif edge[:t_connections] == [2, 3] || edge[:t_connections] == [3, 2]
                            edge[:edge_color] = :orange
                        elseif edge[:t_connections] == [3, 1] || edge[:t_connections] == [1, 3]
                            edge[:edge_color] = :magenta
                        else
                            Pliers.warning_text("Unexpected connections: $(edge[:t_connections])")
                            edge[:edge_color] = :gray
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
                            Pliers.warning_text("Unexpected connections: $(edge[:t_connections])")
                            edge[:edge_color] = :gray
                        end
                    elseif length(edge[:t_connections]) == 4
                        if length(edge[:t_connections]) == 4
                            edge[:edge_color] = :purple
                        else
                            Pliers.warning_text("Unexpected connections: $(edge[:t_connections])")
                            edge[:edge_color] = :gray
                        end
                    else
                        edge[:edge_color] = :gray
                    end
                else
                    edge[:edge_color] = :gray

                end
            end



        end
    else
        # MATH related formatting
        for (_, edge) in network_graph.eprops
            # Set default arrow properties
            edge[:arrow_show] = false
            edge[:arrow_marker] =  '⠀'  # '→' '⠀'
            edge[:arrow_size] = 0
            edge[:arrow_shift] = 0.5

            # Determine if branch works as transformer
            is_transformer = false
            if haskey(edge, :is_transformer) && edge[:is_transformer]
                is_transformer = true
            end

            # Determine if branch is a switch
            is_switch = haskey(edge, :is_switch) && edge[:is_switch]

            # Assume gold for transformers (similar to lines without t_connections in ENG logic)
            if is_transformer
                edge[:edge_color] = :gold
                edge[:arrow_show] = true
                edge[:arrow_marker] = 'Ꝏ'
                edge[:arrow_size] = 20
                edge[:arrow_shift] = 0.5
            elseif is_switch
                # Switch decoration: teal color, marker reflects open/closed state
                edge[:edge_color] = :teal
                edge[:arrow_show] = true
                edge[:arrow_marker] = edge[:state] == 0 ? '□' : '■'
                edge[:arrow_size] = 20
                edge[:arrow_shift] = 0.5
            else
                # Colored depending on phases (size of matrix)
                edge[:edge_color] = :black
                n_phases = 1
                if haskey(edge, :t_connections)
                    n_phases = length(edge[:t_connections])
                end
                if n_phases == 1
                if edge[:t_connections] == [1]
                    edge[:edge_color] = :red
                elseif edge[:t_connections] == [2]
                    edge[:edge_color] = :green
                elseif edge[:t_connections] == [3]
                    edge[:edge_color] = :blue
                elseif edge[:t_connections] == [4]
                    edge[:edge_color] = :black
                else
                    Pliers.warning_text("Unexpected connections: $(edge[:t_connections])")
                    edge[:edge_color] = :gray
                end
            elseif n_phases == 2
                if edge[:t_connections] == [1, 4] || edge[:t_connections] == [4, 1]
                    edge[:edge_color] = :red
                elseif edge[:t_connections] == [2, 4] || edge[:t_connections] == [4, 2]
                    edge[:edge_color] = :green
                elseif edge[:t_connections] == [3, 4] || edge[:t_connections] == [4, 3]
                    edge[:edge_color] = :blue
                elseif edge[:t_connections] == [1, 2] || edge[:t_connections] == [2, 1]
                    edge[:edge_color] = :cyan
                elseif edge[:t_connections] == [2, 3] || edge[:t_connections] == [3, 2]
                    edge[:edge_color] = :orange
                elseif edge[:t_connections] == [3, 1] || edge[:t_connections] == [1, 3]
                    edge[:edge_color] = :magenta
                else
                    Pliers.warning_text("Unexpected connections: $(edge[:t_connections])")
                    edge[:edge_color] = :gray
                end
            elseif n_phases == 3
                if edge[:t_connections] == [1, 2, 3]
                    edge[:edge_color] = :purple
                elseif edge[:t_connections] == [1, 2, 4]
                    edge[:edge_color] = :blue
                elseif edge[:t_connections] == [2, 3, 4]
                    edge[:edge_color] = :blue
                elseif edge[:t_connections] == [3, 1, 4]
                    edge[:edge_color] = :blue
                else
                    Pliers.warning_text("Unexpected connections: $(edge[:t_connections])")
                    edge[:edge_color] = :gray
                end
            elseif n_phases == 4
                if length(edge[:t_connections]) == 4
                    edge[:edge_color] = :purple
                else
                    Pliers.warning_text("Unexpected connections: $(edge[:t_connections])")
                    edge[:edge_color] = :gray
                end
            else
                @warn "Unexpected number of phases: $n_phases for edge with connections $(get(edge, :t_connections, "N/A")), defaulting to gray"
                edge[:edge_color] = :gray

            end
            end

            
        end
    end
end

# function _decorate_edges!(network_graph::MetaDiGraph, data::Dict{String,Any}, phase::String)
#     if _is_eng(data)
#         for (_, edge) in network_graph.eprops
#             if !haskey(edge, :t_connections) # I am suing this to determine that it is a transformer not a line
#                 edge[:edge_color] = :gold
#             else
#                 if !isempty(edge[:t_connections])
#                     if length(edge[:t_connections]) == 1
#                         if edge[:t_connections] == [1]
#                             edge[:edge_color] = :red
#                         elseif edge[:t_connections] == [2]
#                             edge[:edge_color] = :white
#                         elseif edge[:t_connections] == [3]
#                             edge[:edge_color] = :white
#                         elseif edge[:t_connections] == [4]
#                             edge[:edge_color] = :white
#                         else
#                             error("Unexpected connections: $(edge[:t_connections])")
#                         end
#                     elseif length(edge[:t_connections]) == 2
#                         if edge[:t_connections] == [1, 4]
#                             edge[:edge_color] = :red
#                         elseif edge[:t_connections] == [2, 4]
#                             edge[:edge_color] = :white
#                         elseif edge[:t_connections] == [3, 4]
#                             edge[:edge_color] = :white
#                         elseif edge[:t_connections] == [1, 2]
#                             edge[:edge_color] = :white
#                         elseif edge[:t_connections] == [2, 3]
#                             edge[:edge_color] = :white
#                         elseif edge[:t_connections] == [3, 1]
#                             edge[:edge_color] = :white
#                         else
#                             error("Unexpected connections: $(edge[:t_connections])")
#                         end
#                     elseif length(edge[:t_connections]) == 3
#                         if edge[:t_connections] == [1, 2, 3]
#                             edge[:edge_color] = :red
#                         elseif edge[:t_connections] == [1, 2, 4]
#                             edge[:edge_color] = :red
#                         elseif edge[:t_connections] == [2, 3, 4]
#                             edge[:edge_color] = :white
#                         elseif edge[:t_connections] == [3, 1, 4]
#                             edge[:edge_color] = :red
#                         else
#                             error("Unexpected connections: $(edge[:t_connections])")
#                         end
#                     elseif length(edge[:t_connections]) == 4
#                         if edge[:t_connections] == [1, 2, 3, 4]
#                             edge[:edge_color] = :red
#                         else
#                             error("Unexpected connections: $(edge[:t_connections])")
#                         end
#                     else
#                         edge[:edge_color] = :white
#                     end
#                 else
#                     edge[:edge_color] = :white

#                 end
#             end
#         end
#     else
#         #TODO: MATH related formatting
#     end
# end

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
    return [string(props(network_graph, i)[:bus_id]) for i in 1:nv(network_graph)]
end

## LOADS
function _write_lolabels(network_graph, data)
    lolabels = Vector{String}(undef, nv(network_graph))
    for i in 1:nv(network_graph)
        loads = props(network_graph, i)[:loads]
        if !isempty(loads)
            lolabels[i] = join([string(load[:load_id]) for load in loads], ", ")
        else
            lolabels[i] = ""
        end
    end
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









# ---------------------------------------------------------------------------
# Voltage-level plotting helpers
# ---------------------------------------------------------------------------

"""
    _extract_vbase_map_eng(eng_sym::Dict) -> Dict{Symbol, Float64}

Build a bus → base-voltage (kV) mapping from an engineering model.

**Strategy**
1. Read `vm_nom` from every transformer winding → direct bus-to-kV assignment.
2. Read `vm_nom` from voltage sources (if present).
3. Propagate known kV values through line and switch connections so that all
   buses in the same voltage zone are labelled, not just transformer terminals.

This covers the common case where engineering-model bus dicts contain no
`vbase` or `vm_nom` key of their own.
"""
function _extract_vbase_map_eng(eng_sym::Dict)
    vbase_map = Dict{Symbol,Float64}()

    @debug "_extract_vbase_map_eng: building vbase map from engineering model"

    # --- Step 1: Transformer windings ---
    if haskey(eng_sym, :transformer)
        @debug "  Found $(length(eng_sym[:transformer])) transformer(s)"
        for (tid, trans) in eng_sym[:transformer]
            @debug "  Transformer '$tid' keys: $(sort(string.(collect(keys(trans)))))"
            buses   = get(trans, :bus,    String[])
            vm_nom  = get(trans, :vm_nom, Float64[])
            @debug "    buses=$buses  vm_nom=$vm_nom"
            for (i, bus) in enumerate(buses)
                if i <= length(vm_nom) && !ismissing(vm_nom[i])
                    vbase_map[Symbol(bus)] = Float64(vm_nom[i])
                    @debug "    → bus :$(Symbol(bus)) assigned $(vm_nom[i]) kV"
                end
            end
        end
    else
        @debug "  No transformers found in engineering model"
    end

    # --- Step 2: Voltage sources ---
    if haskey(eng_sym, :voltage_source)
        @debug "  Found $(length(eng_sym[:voltage_source])) voltage source(s)"
        for (vsid, vs) in eng_sym[:voltage_source]
            @debug "  VoltageSource '$vsid' keys: $(sort(string.(collect(keys(vs)))))"
            bus_id = Symbol(get(vs, :bus, ""))
            if haskey(vs, :vm_nom) && !ismissing(vs[:vm_nom])
                vn = vs[:vm_nom]
                val = Float64(isa(vn, AbstractArray) ? first(vn) : vn)
                vbase_map[bus_id] = val
                @debug "    → source bus :$bus_id assigned $val kV from vm_nom"
            else
                @debug "    VoltageSource '$vsid' has no vm_nom key"
            end
        end
    end

    # --- Step 3: Propagate through lines (same voltage zone) ---
    for component_key in (:line, :switch)
        if haskey(eng_sym, component_key)
            @debug "  Propagating through $(length(eng_sym[component_key])) $(component_key)(s)"
            changed = true
            pass = 0
            while changed
                changed = false
                pass += 1
                for (_, elem) in eng_sym[component_key]
                    f = Symbol(get(elem, :f_bus, ""))
                    t = Symbol(get(elem, :t_bus, ""))
                    if haskey(vbase_map, f) && !haskey(vbase_map, t)
                        vbase_map[t] = vbase_map[f]
                        @debug "    pass $pass: :$t ← :$f ($(vbase_map[f]) kV)"
                        changed = true
                    elseif haskey(vbase_map, t) && !haskey(vbase_map, f)
                        vbase_map[f] = vbase_map[t]
                        @debug "    pass $pass: :$f ← :$t ($(vbase_map[t]) kV)"
                        changed = true
                    end
                end
            end
        end
    end

    @debug "_extract_vbase_map_eng result: $vbase_map"
    return vbase_map
end


"""
    _extract_vbase_map_math(math_sym::Dict) -> Dict{Symbol, Float64}

Build a bus → base-voltage (kV) mapping from a mathematical model.
Math-model buses carry a `vbase` field populated by `transform_data_model`.
"""
function _extract_vbase_map_math(math_sym::Dict)
    vbase_map = Dict{Symbol,Float64}()

    @debug "_extract_vbase_map_math: reading vbase from math bus dicts"

    if !haskey(math_sym, :bus)
        @debug "  No :bus key in math model"
        return vbase_map
    end

    for (bid, bus) in math_sym[:bus]
        @debug "  Bus '$bid' keys: $(sort(string.(collect(keys(bus)))))"
        if haskey(bus, :vbase) && !ismissing(bus[:vbase])
            vbase_map[Symbol(bid)] = Float64(bus[:vbase])
            @debug "    → bus :$bid vbase=$(bus[:vbase]) kV"
        else
            @debug "    → bus :$bid has no vbase key"
        end
    end

    @debug "_extract_vbase_map_math result: $vbase_map"
    return vbase_map
end


"""
    _get_bus_vbase(bus_id::Symbol, vbase_map::Dict) -> Union{Float64, Symbol}

Look up the base voltage for `bus_id` in `vbase_map`.
Returns `:Unknown` when the bus is not in the map.
"""
function _get_bus_vbase(bus_id::Symbol, vbase_map::Dict)
    if haskey(vbase_map, bus_id) && !ismissing(vbase_map[bus_id])
        v = vbase_map[bus_id]
        @debug "  _get_bus_vbase(:$bus_id) → $v kV"
        return Float64(v)
    end
    @debug "  _get_bus_vbase(:$bus_id) → :Unknown (not in vbase_map)"
    return :Unknown
end


"""
    _format_voltage_label(v) -> String

Return a human-readable voltage label.
Values ≥ 1 are assumed to be in kV; smaller values are converted to V.
`:Unknown` is returned as `"Unknown"`.
"""
function _format_voltage_label(v)
    v === :Unknown && return "Unknown"
    v_kv = Float64(v)
    if v_kv >= 1.0
        return "$(round(v_kv, digits=3)) kV"
    else
        return "$(round(v_kv * 1000.0, digits=1)) V"
    end
end


"""
    _build_voltage_colormap(voltage_values) -> Dict

Build a mapping from unique voltage values to distinct Makie colours.

Known voltages receive colours from the Wong colour-blind-safe palette,
cycling if there are more levels than palette entries.  `:Unknown` maps
to `:gray`.
"""
function _build_voltage_colormap(voltage_values)
    @debug "_build_voltage_colormap: input values = $voltage_values"
    known = sort(unique([v for v in voltage_values if v !== :Unknown]))
    @debug "  distinct known voltages (sorted): $known"
    palette = Makie.wong_colors()
    color_map = Dict{Any,Any}()
    for (i, v) in enumerate(known)
        color_map[v] = palette[mod1(i, length(palette))]
        @debug "  color_map[$v] = $(color_map[v])"
    end
    if :Unknown in voltage_values
        color_map[:Unknown] = :gray
        @debug "  color_map[:Unknown] = :gray"
    end
    @debug "_build_voltage_colormap result: $(keys(color_map))"
    return color_map
end


"""
    _decorate_nodes_by_voltage!(network_graph, data, vbase_map) -> (Dict, Vector, Vector, Vector)

Set `:node_color`, `:node_marker`, and `:marker_size` on every vertex of
`network_graph` according to the bus base voltage looked up from `vbase_map`.

Returns `(voltage_color_map, node_color, node_marker, node_size)` where
`voltage_color_map` maps each voltage value to its assigned colour.
"""
function _decorate_nodes_by_voltage!(
    network_graph::MetaDiGraph,
    data::Dict{String,Any},
    vbase_map::Dict,
)
    @debug "_decorate_nodes_by_voltage!: decorating $(nv(network_graph)) nodes"

    # Collect per-node vbase values keyed by vertex index for fast lookup
    node_vbase = Dict{Int,Any}()
    for i in 1:nv(network_graph)
        node = props(network_graph, i)
        bus_id = get(node, :bus_id, Symbol(""))
        vb = _get_bus_vbase(bus_id, vbase_map)
        node_vbase[i] = vb
        @debug "  vertex $i (:$bus_id) → vbase = $vb"
    end

    all_vbase = collect(values(node_vbase))
    @debug "  all collected vbase values: $all_vbase"
    color_map = _build_voltage_colormap(all_vbase)

    for i in 1:nv(network_graph)
        node = props(network_graph, i)
        vb = node_vbase[i]
        node[:node_color]  = color_map[vb]
        node[:node_marker] = :circle
        node[:marker_size] = 10
        @debug "  vertex $i color=$(color_map[vb])  (vbase=$vb)"
    end

    node_color  = [props(network_graph, i)[:node_color]  for i in 1:nv(network_graph)]
    node_marker = [props(network_graph, i)[:node_marker] for i in 1:nv(network_graph)]
    node_size   = [props(network_graph, i)[:marker_size]  for i in 1:nv(network_graph)]

    @debug "_decorate_nodes_by_voltage!: done. color_map keys = $(keys(color_map))"
    return color_map, node_color, node_marker, node_size, node_vbase
end


"""
    _decorate_edges_by_voltage!(network_graph, node_color, node_vbase, color_map)

Colour each edge in `network_graph` using the voltage zone of its **source**
vertex (`src(e)`).  Edges that cross voltage zones (transformers) use the
source-side colour.  All arrow decorations are suppressed.
"""
function _decorate_edges_by_voltage!(
    network_graph::MetaDiGraph,
    node_color::Vector,
    node_vbase::Dict,
    color_map::Dict,
)
    @debug "_decorate_edges_by_voltage!: decorating $(ne(network_graph)) edges"

    for e in edges(network_graph)
        s = src(e)
        d = dst(e)
        src_vb  = get(node_vbase, s, :Unknown)
        dst_vb  = get(node_vbase, d, :Unknown)
        src_col = node_color[s]

        if src_vb != dst_vb
            @debug "  edge ($s→$d): crosses voltage zones ($src_vb kV → $dst_vb kV), using src colour"
        else
            @debug "  edge ($s→$d): same zone ($src_vb kV)"
        end

        set_prop!(network_graph, e, :edge_color,   src_col)
        set_prop!(network_graph, e, :arrow_show,   false)
        set_prop!(network_graph, e, :arrow_marker, '⠀')
        set_prop!(network_graph, e, :arrow_size,   0)
        set_prop!(network_graph, e, :arrow_shift,  0.5)
    end
end


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

"""
    plot_network_by_voltage(data::Dict{String,Any}; kwargs...) -> Figure

Plot a network graph where **node and edge colours reflect the base voltage
level** (`v_base`) of each bus, rather than the phase-load decoration used
by [`plot_network_tree`](@ref).

A legend mapping each colour to its voltage level (e.g. "0.4 kV", "11 kV")
is appended to the right of the plot.  Buses without a detectable `vbase`
are shown in gray and labelled *"Unknown"*.

Both engineering model (`data_model = ENGINEERING`) and mathematical model
(`data_model = MATHEMATICAL`) are supported.  For the engineering model the
voltage base is inferred from transformer `vm_nom` fields and propagated
through line connections.  For the mathematical model it is read directly
from the `vbase` field of each bus.

Enables `@debug` logging (`ENV["JULIA_DEBUG"] = "Pliers"` or
`Logging.with_logger(Logging.ConsoleLogger(stderr, Logging.Debug))`) to
inspect per-bus vbase resolution.

# Arguments
- `data::Dict{String,Any}`: Engineering or Mathematical network model.
- `figure_size`: Overall figure resolution in pixels. Default `(1200, 900)`.
- `show_node_labels`: Show bus-id labels on nodes. Default `false`.
- `show_edge_labels`: Show line-id labels on edges. Default `false`.
- `show_load_labels`: Show load labels on nodes. Default `false`.
- `layout`: Node layout algorithm. Default [`smart_layout`](@ref).
- `makie_backend`: Makie backend module. Default `WGLMakie`.
- `kwargs...`: Additional keyword arguments forwarded to `graphplot!`.

# Returns
A `Makie.Figure` containing the network plot and the voltage legend.

# Example
```julia
using WGLMakie
using Logging
# optional: enable debug output
Logging.global_logger(Logging.ConsoleLogger(stderr, Logging.Debug))

eng = PowerModelsDistribution.parse_file("network.dss")
fig = plot_network_by_voltage(eng)
```
"""
function plot_network_by_voltage(
    data::Dict{String,Any};
    figure_size=(1200, 900),
    show_node_labels=false,
    show_edge_labels=false,
    show_load_labels=false,
    layout=smart_layout,
    edge_labels_type=:line_id,
    phase="1",
    makie_backend=WGLMakie,
    kwargs...
)
    makie_backend.activate!()

    @debug """plot_network_by_voltage: data_model=$(get(data, "data_model", "(unknown)"))"""
    # --- Build the network meta-graph ---
    network_graph, _, _ = create_network_graph(data, layout)
    @debug "  graph: $(nv(network_graph)) vertices, $(ne(network_graph)) edges"

    # --- Build vbase map depending on model type ---
    data_sym = convert_keys_to_symbols(data)
    if _is_eng(data)
        @debug "  Using engineering-model vbase extraction"
        vbase_map = _extract_vbase_map_eng(data_sym)
    else
        @debug "  Using mathematical-model vbase extraction"
        vbase_map = _extract_vbase_map_math(data_sym)
    end

    @debug "  vbase_map ($(length(vbase_map)) entries): $vbase_map"

    # --- Labels ---
    nlabels = show_node_labels ? _write_nlabels(network_graph, data) : nothing
    lolabels = show_load_labels ? _write_lolabels(network_graph, data) : nothing
    elabels  = show_edge_labels ?
        (edge_labels_type == :line_id ?
            _write_line_id_elabels(network_graph, data) :
            _write_results_elabels(network_graph, data, phase)) :
        nothing

    if show_load_labels
        if show_node_labels
            nlabels = [string(n, "\n", l) for (n, l) in zip(nlabels, lolabels)]
        else
            nlabels = lolabels
        end
    end

    # --- Node decoration: voltage-based colours ---
    color_map, node_color, node_marker, node_size, node_vbase =
        _decorate_nodes_by_voltage!(network_graph, data, vbase_map)

    # --- Edge decoration: voltage-zone colours (matches node colour) ---
    _decorate_edges_by_voltage!(network_graph, node_color, node_vbase, color_map)
    edge_color   = [get_prop(network_graph, e, :edge_color)   for e in edges(network_graph)]
    arrow_marker = [get_prop(network_graph, e, :arrow_marker) for e in edges(network_graph)]
    arrow_size   = [get_prop(network_graph, e, :arrow_size)   for e in edges(network_graph)]
    arrow_shift  = [get_prop(network_graph, e, :arrow_shift)  for e in edges(network_graph)]

    # --- Build figure with legend column ---
    fig = Figure(size=figure_size)
    ax  = Axis(fig[1, 1])

    graphplot!(
        ax,
        network_graph;
        layout=layout,
        nlabels=nlabels,
        show_node_labels=show_node_labels || show_load_labels,
        elabels=elabels,
        show_edge_labels=show_edge_labels,
        node_color=node_color,
        node_marker=node_marker,
        node_size=node_size,
        edge_color=edge_color,
        arrow_show=true,
        arrow_marker=arrow_marker,
        arrow_size=arrow_size,
        arrow_shift=arrow_shift,
        edge_plottype=:linesegments,
        kwargs...
    )

    # --- Legend: sorted known voltages first, Unknown last ---
    known_voltages  = sort([v for v in keys(color_map) if v !== :Unknown])
    has_unknown     = :Unknown in keys(color_map)
    sorted_voltages = vcat(known_voltages, has_unknown ? [:Unknown] : [])

    legend_elements = [
        MarkerElement(color=color_map[v], marker=:circle, markersize=15)
        for v in sorted_voltages
    ]
    legend_labels = [_format_voltage_label(v) for v in sorted_voltages]

    @debug "  legend entries: $(collect(zip(legend_labels, sorted_voltages)))"

    Legend(
        fig[1, 2],
        legend_elements,
        legend_labels,
        "Voltage Level";
        framevisible=true,
        tellheight=false,
    )

    colsize!(fig.layout, 1, Relative(0.82))

    return fig
end


# Re-export plotting functions from parent module
export create_graph
export plot_network_tree, plot_network_tree!
export plot_network_coords, plot_network_coords!
export plot_network_map
export plot_network_by_voltage
export bus_phasor, bus_phasor!
export smart_layout
end # module PMDGraph
