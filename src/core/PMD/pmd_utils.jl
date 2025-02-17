
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




### Operation on Results Dictionary ###
"""
    fluff_bus_voltages(PF_Res::Dict{String, Any})

This function is used to get all possible voltage forms from the output of the OPF voltage variables, so it processes the bus voltages in the given Power Flow Results dictionary `PF_Res`. For each bus in `PF_Res["bus"]`, it checks if the bus has real (`vr`) and imaginary (`vi`) voltage components. If both components are present, it calculates the complex voltage `V`, the voltage magnitude `vm`, and the voltage angle `va` for the bus. A warning is issued indicating that only `vr` and `vi` can be processed.

# Arguments
- `PF_Res::Dict{String, Any}`: A dictionary containing the Power Flow Results, including bus voltage information.

# Modifies
- Adds the following keys to each bus dictionary if `vr` and `vi` are present:
  - `V`: The complex voltage calculated as `vr + vi*im`.
  - `vm`: The magnitude of the complex voltage.
  - `va`: The angle of the complex voltage.

"""
function fluff_bus_voltages!(PF_Res::Dict{String, Any})
    
    for (_, bus) in PF_Res["solution"]["bus"]
        if  haskey(bus, "vr") && haskey(bus, "vi") 
            bus["V"] = bus["vr"] .+ bus["vi"]*im
            bus["vm"] = abs.(bus["V"])
            bus["va"] = angle.(bus["V"])
        elseif haskey(bus, "va") && haskey(bus, "vm")
            bus["V"] = bus["vm"] .* exp.(im.*bus["va"])
            bus["vr"] = real.(bus["V"])
            bus["vi"] = imag.(bus["V"])
        else 
            warning_text("Can't fluff $(keys(bus)). Now I can just fluff the (`vr`, `vi`) and (`va`, `vm`) pairs.")
        end 
    end

end 


function dictify_voltages!(PF_res::Dict{String, Any}, math)
    
    fluff_bus_voltages!(PF_res)
    pf_sol = PF_res["solution"]

    for (b, bus) in pf_sol["bus"]
        terminals = math["bus"][b]["terminals"]
        # write a dictionary where the key is the terminal number and the value is the voltage at that terminal
        bus["term"] = Dict(string(term) => bus["V"][i] for (i, term) in enumerate(terminals))
    end 

end