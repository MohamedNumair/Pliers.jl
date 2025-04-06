
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

    if haskey(pf_sol["bus"][bus_index]["voltage"], string(_N_IDX))  
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
        math_meas["bus"][b]["terminals"] = setdiff(math["bus"][b]["terminals"], _N_IDX)
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

