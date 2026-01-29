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
network_graph, layout, eng_sym = create_network_graph_eng(eng, GraphMakie.Buchheim())
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
        line[:linecodes] = eng_sym[:linecode][Symbol(line[:linecode])]

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
network_graph, layout, math_sym = create_network_graph_math(math, GraphMakie.Buchheim())
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