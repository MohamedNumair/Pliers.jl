
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

    for (b, bus) in pf_sol["bus"]
        terminals = math["bus"][b]["terminals"]
        # write a dictionary where the key is the terminal number and the value is the voltage at that terminal
        bus["voltage"] = Dict(string(term) => bus["V"][i] for (i, term) in enumerate(terminals))
    end 

end

function solution_dictify_loads!(pf_sol::Dict{String, Any}, math::Dict{String, Any}; formulation = "IVR") 
    # The idea is to create the complex current and power for each load and store it in the load dictionary under the key "current" and "power" respectively.
    # The current is calculated as `I = P + Q*im` and the power is calculated as `S = P + Q*im`  
    for (l, load) in pf_sol["load"]
        terminals = math["load"][l]["connections"]
        # write a dictionary where the key is the terminal number and the value is the current at that terminal
        load["current"] = haskey(load,"crd_bus") ? Dict(string(term) => load["crd_bus"][i] + load["cid_bus"][i]*im for (i, term) in enumerate(terminals)) : Dict(string(term) => load["crd"][i] + load["cid"][i]*im for (i, term) in enumerate(terminals)) 
        haskey(load, "pd") || haskey(load, "pd_bus") ? load["power"] =  haskey(load,"pd_bus") ?  Dict(string(term) => load["pd_bus"][i] + load["qd_bus"][i]*im for (i, term) in enumerate(terminals)) :  Dict(string(term) => load["pd"][i] + load["qd"][i]*im for (i, term) in enumerate(terminals))  : nothing
    end
end 

function solution_dictify_branches!(pf_sol::Dict{String, Any}, math::Dict{String, Any}; formulation = "IVR")
    # The idea is to create the complex current and power for each branch and store it in the branch dictionary under the key "current" and "power" respectively.
    # The current is calculated as `I = P + Q*im` and the power is calculated as `S = P + Q*im`  
    for (b, branch) in pf_sol["branch"]
        f_terminals = math["branch"][b]["f_connections"]
        t_terminals = math["branch"][b]["t_connections"]
        # write a dictionary where the key is the terminal number and the value is the current at that terminal
        branch["current_from"] = Dict(string(term) => branch["cr_fr"][i] + branch["ci_fr"][i]*im for (i, term) in enumerate(f_terminals))
        haskey(branch, "csr_fr") ? branch["shunt_current_from"] = Dict(string(term) => branch["csr_fr"][i] + branch["csi_fr"][i]*im for (i, term) in enumerate(f_terminals)) : nothing
        branch["current_to"] = Dict(string(term) => branch["cr_to"][i] + branch["ci_to"][i]*im for (i, term) in enumerate(t_terminals))
        haskey(branch, "csr_to") ? branch["shunt_current_to"] = Dict(string(term) => branch["csr_to"][i] + branch["csi_to"][i]*im for (i, term) in enumerate(t_terminals)) : nothing
        branch["power_from"] = Dict(string(term) => branch["pf"][i] + branch["qf"][i]*im for (i, term) in enumerate(f_terminals))
        branch["power_to"] = Dict(string(term) => branch["pt"][i] + branch["qt"][i]*im for (i, term) in enumerate(t_terminals))
    end
end

function solution_dictify_gens!(pf_sol::Dict{String, Any}, math::Dict{String, Any}; formulation = "IVR")
    # The idea is to create the complex current and power for each generator and store it in the generator dictionary under the key "current" and "power" respectively.
    # The current is calculated as `I = P + Q*im` and the power is calculated as `S = P + Q*im`  
    for (g, gen) in pf_sol["gen"]
        terminals = setdiff(math["gen"][g]["connections"], [4]) 
        # write a dictionary where the key is the terminal number and the value is the current at that terminal
        gen["current"] = Dict(string(term) => gen["crg"][i] + gen["cig"][i]*im for (i, term) in enumerate(terminals))
        haskey(gen, "pg") || haskey(gen, "pg_bus") ? gen["power"] = Dict(string(term) => gen["pg"][i] + gen["qg"][i]*im for (i, term) in enumerate(terminals)) : nothing
    end
end

function dictify_solution!(pf_sol::Dict{String, Any}, math::Dict{String, Any}; formulation = "IVR")
    solution_dictify_buses!(pf_sol, math; formulation = formulation)
    solution_dictify_loads!(pf_sol, math; formulation = formulation)
    solution_dictify_branches!(pf_sol, math; formulation = formulation)
    solution_dictify_gens!(pf_sol, math; formulation = formulation)
end

"""
    separate_phase_neutral_voltages(pf_sol, bus_index)

Separate the phase and neutral voltages from the power flow solution for a given bus.

# Arguments
- `pf_sol::Dict`: The power flow solution dictionary containing voltage information.
- `bus_index::Int`: The index of the bus for which to separate the voltages.

# Returns
- `phase_voltage::Vector{ComplexF64}`: A vector containing the phase voltages (up to 3 phases).
- `neutral_voltage::ComplexF64`: The neutral voltage. If the neutral voltage is not present, returns 0 + 0im "assuming grounded".

# Notes
- The neutral terminal is currently hardcoded to "4". This should be fixed in future versions.
"""

function separate_phase_neutral_voltages(pf_sol, bus_index)
    phase_voltage = []
    for i in 1:3 
        if haskey(pf_sol["bus"][bus_index]["voltage"], string(i))
            push!(phase_voltage, pf_sol["bus"][bus_index]["voltage"][string(i)])
        end
    end

    neutral_voltage = ComplexF64[]

    if haskey(pf_sol["bus"][bus_index]["voltage"], "4")  #TODO: fix hardcoding the neutral terminal
        neutral_voltage = pf_sol["bus"][bus_index]["voltage"]["4"]
    else 
        neutral_voltage = 0 + 0im # assuming grounded
    end

    return phase_voltage, neutral_voltage
end


function phase_neutral_voltage_file(math, pf_sol)
    math_meas = deepcopy(math)
    for (b, bus) in pf_sol["bus"]
        pvs, vn = Pliers.separate_phase_neutral_voltages(pf_sol, b)
        vmn = abs.(pvs .- vn)
        bus["vmn"] = vmn
        math_meas["bus"][b]["terminals"] = setdiff(math["bus"][b]["terminals"], 4)
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

function write_delta_readings(math, pf_sol)
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

