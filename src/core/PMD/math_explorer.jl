
# MATH report d=====(￣▽￣*)b

function math_report(math::Dict{String, Any}; detailed = false)
    print(UNDERLINE(BLUE_FG("Report for the ",   BOLD("$(get(math, "data_model", nothing))"),
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
        index = Int[],
        bus_id = String[],
        bus_type = Int[],
        terminals = Vector{Int}[], 
        grounded = Vector{Int}[],
        vm_pair_ub = Vector{Tuple{Any, Any, Real}}[],
        vm_pair_lb = Vector{Tuple{Any, Any, Real}}[],
        source_id = String[],
        vbase = Float64[],
        vmin = Vector{Float64}[],
        vmax = Vector{Float64}[],
     )
end

function _pushing_math_buses_df!(buses_df, bus_id, bus)
    push!(buses_df, (   
        get(bus, "index", missing),
        bus_id,
        get(bus, "bus_type", missing),
        get(bus, "terminals", Int[]),
        get(bus, "grounded", Int[]),
        get(bus, "vm_pair_ub", Tuple{Any, Any, Real}[]),
        get(bus, "vm_pair_lb", Tuple{Any, Any, Real}[]),
        get(bus, "source_id", missing),
        get(bus, "vbase", missing),
        get(bus, "vmin", Float64[]),
        get(bus, "vmax", Float64[]),
    ))
end

function math_buses_table(math::Dict{String, Any})
    buses = haskey(math, "bus") ? math["bus"] : error("No buses found in the MATHEMATICAL data check the data model has \"bus\" key") 
    buses_df = _initiate_math_buses_df()  
    for (bus_id, bus) in buses
        _pushing_math_buses_df!(buses_df, bus_id, bus)
    end
    header("Buses Table ($(nrow(buses_df)) buses)")
    pretty_table(sort!(buses_df))
    extra_keys(buses, ["index", "bus_id", "bus_type", "terminals", "grounded", "vm_pair_ub", "vm_pair_lb", "source_id", "vbase", "vmin", "vmax"])
end

function math_buses_table(math::Dict{String, Any}, condition::Function)
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
        table = DataFrame(labels = labels, values = values)
        pretty_table(table, header=["Name", "Value"])
    end
end

# MATH BRANCHES (. ❛ ᴗ ❛.) 

function _intitiate_math_branches_df()
    return DataFrame(   
                        index = Int[],
                        source_id = String[],
                        f_bus = Int[],
                        t_bus = Int[],
                        f_connections = Vector{Int}[],
                        t_connections = Vector{Int}[],
                        br_status = Int[],
                        vbase = Float64[],
                        br_r = Matrix{Float64}[],
                        br_x = Matrix{Float64}[],
                        g_fr = Matrix{Float64}[],
                        b_fr = Matrix{Float64}[],
                        g_to = Matrix{Float64}[],
                        b_to = Matrix{Float64}[],
                        rate_a = Vector{Real}[],
                        rate_b = Vector{Real}[],
                        rate_c = Vector{Real}[],
                        c_rating_a = Vector{Real}[],
                        c_rating_b = Vector{Real}[],
                        c_rating_c = Vector{Real}[],
                        angmin = Vector{Real}[],
                        angmax = Vector{Real}[],
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


function math_branches_table(math::Dict{String, Any})
    branches = haskey(math, "branch") ? math["branch"] : error("No branches found in the MATHEMATICAL data check the data model has \"branch\" key") 
    branches_df = _intitiate_math_branches_df()
    for (branch_id, branch) in branches
        _pushing_math_branches_df!(branches_df, branch_id, branch)
    end
    header("Branches Table ($(nrow(branches_df)) branches)")
    pretty_table(sort!(branches_df))
    extra_keys(branches, ["index", "name", "source_id", "f_bus", "t_bus", "f_connections", "t_connections", "br_status", "vbase", "br_r", "br_x", "g_fr", "b_fr", "g_to", "b_to", "rate_a", "rate_b", "rate_c", "c_rating_a", "c_rating_b", "c_rating_c", "angmin", "angmax"])
end 

function math_branches_table(math::Dict{String, Any}, condition::Function)
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

function math_branch_details(math::Dict{String, Any}, branch_idx)
    branches = haskey(math, "branch") ? math["branch"] : error("No branches found in the MATHEMATICAL data check the data model has \"branch\" key") 
    for id in branch_idx
        branch = branches[id]
        
        header("Branch $id")

        labels_1 = ["f_bus", "t_bus", "f_connections", "t_connections", "br_status", "vbase","rate_a", "rate_b", "rate_c", "c_rating_a", "c_rating_b", "c_rating_c", "angmin", "angmax"]
        values_1 = [get(branch, label, missing) for label in labels_1]
        table_1 = DataFrame(labels = labels_1, values = values_1)
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
        index = Int[],
        source_id = String[],
        load_bus = Int[],
        connections = Vector{Int}[], 
        configuration = String[],
        name = String[],
        pd = Vector{Real}[],
        qd = Vector{Real}[],
        status = Int[],
        vbase = Float64[],
        vnom_kv = Float64[],
        dispatchable = Int[],
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

function math_loads_table(math::Dict{String, Any})
    loads = haskey(math, "load") ? math["load"] : error("No loads found in the MATHEMATICAL data check the data model has \"load\" key") 
    loads_df = _initiate_math_loads_df()
    for (load_id, load) in loads
        _pushing_math_loads_df!(loads_df, load_id, load)
    end
    header("Loads Table ($(nrow(loads_df)) loads)")
    pretty_table(sort!(loads_df))
    extra_keys(loads, ["index", "source_id", "load_bus", "connections", "configuration", "name", "status", "qd", "vbase", "vnom_kv", "dispatchable", "pd"])       
end
function math_loads_table(math::Dict{String, Any}, condition::Function)
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
        table = DataFrame(labels = labels, values = values)
        pretty_table(table, header=["Name", "Value"])
    end
end


# MATH GEN (. ❛ ᴗ ❛.)


function _initiate_math_gen_df()
    return DataFrame(   
                        index = Int[],
                        gen_bus = Int[],
                        source_id = String[],
                        connections = Vector{Int}[],
                        configuration = String[],
                        name = String[],
                        pg = Vector{Real}[],
                        qg = Vector{Real}[],
                        pmax = Vector{Real}[],
                        pmin = Vector{Real}[],
                        qmax = Vector{Real}[],
                        qmin = Vector{Real}[],
                        vg = Vector{Real}[],
                        vbase = Float64[],
                        gen_status = Int[],
                        model = Int[],
                        cost = Vector{Real}[],
                        ncost = Int[],
                        control_mode = Int[],
                        startup = Float64[],
                        shutdown = Float64[],
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

function math_gen_table(math::Dict{String, Any})
    gens = haskey(math, "gen") ? math["gen"] : error("No generators found in the MATHEMATICAL data check the data model has \"gen\" key") 
    gens_df = _initiate_math_gen_df()
    for (gen_id, gen) in gens
        _pushing_math_gen_df!(gens_df, gen_id, gen)
    end
    header("Generators Table ($(nrow(gens_df)) generators)")
    pretty_table(sort!(gens_df))
    extra_keys(gens, ["index", "source_id", "gen_bus", "connections", "configuration", "name", "pg", "qg", "pmax", "pmin", "qmax", "qmin", "vg", "vbase", "gen_status", "model", "cost", "ncost", "control_mode", "startup", "shutdown"])  
end

function math_gen_table(math::Dict{String, Any}, condition::Function)
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
        table = DataFrame(labels = labels, values = values)
        pretty_table(table, header=["Name", "Value"])
    end
end


# MATH TRANSFORMERS (⚡)

function _initiate_math_transformers_df()
    return DataFrame(   
        index = Int[],
        name = String[],
        source_id = String[],
        f_bus = Int[],
        t_bus = Int[],
        f_connections = Vector{Int}[],
        t_connections = Vector{Int}[],
        configuration = String[],
        tm_nom = Float64[],
        tm_set = Vector{Float64}[],
        tm_fix = Vector{Bool}[],
        tm_step = Vector{Float64}[],
        f_vbase = Float64[],
        t_vbase = Float64[],
        polarity = Int[],
        sm_ub = Float64[],
        cm_ub = Float64[],
        status = Int[],
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
function math_transformers_table(math::Dict{String, Any})
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
function math_transformers_table(math::Dict{String, Any}, condition::Function)
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
        table = DataFrame(labels = labels, values = values)
        pretty_table(table, header=["Name", "Value"])
    end
end


