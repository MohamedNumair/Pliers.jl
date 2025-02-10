
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
        if data["data_model"] == PowerModelsDistribution.ENGINEERING
            return true
        elseif data["data_model"] == PowerModelsDistribution.MATHEMATICAL
            return false
        else
            error("Invalid data model found in the provided model dictionary. it has to be either ENGINEERING or MATHEMATICAL")
        end
    else
        error("No data model found in the provided model dictionary.")
    end
end

function _is_eng(graph::MetaDiGraph)

    if haskey(first(graph.eprops).second, :branch_id)
        return false
    else
        return true
    end
    error("I don't know if this graph is for a MATHEMATICAL or ENGINEERING model. I usually check if it has `:branch_id` or `:line_id` in the edge properties to tell which one it is.")
end