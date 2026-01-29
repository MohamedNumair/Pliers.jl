
# API Reference {#API-Reference}

This page lists all exported functions and types from Pliers.jl. 

All functions are available directly from the main `Pliers` module. For convenience, related functions are also re-exported through sub-modules that group them by domain.

## Main Module {#Main-Module}
<details class='jldocstring custom-block' open>
<summary><a id='Pliers.Pliers' href='#Pliers.Pliers'><span class="jlbinding">Pliers.Pliers</span></a> <Badge type="info" class="jlObjectType jlModule" text="Module" /></summary>



```julia
Pliers
```


A Julia package providing tools for analyzing power distribution systems. Designed to be used in conjunction with PowerModelsDistribution.jl and PowerModelsDistributionStateEstimation.jl for simplified reporting, analysis, and visualization.

**Sub-modules**

The package provides optional sub-module access for organized imports:
- `PMDUtils`: Re-exports PMD-related utility functions
  
- `PMDSEUtils`: Re-exports PMDSE-related utility functions  
  
- `PMDPlotting`: Re-exports plotting functions
  

**Author**

Mohamed Numair (mnumair.com)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers._calculate_MAPE-Tuple{Any, Any, Any}' href='#Pliers._calculate_MAPE-Tuple{Any, Any, Any}'><span class="jlbinding">Pliers._calculate_MAPE</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_calculate_MAPE(SE_RES, PF_RES, math)
```


Calculate the Mean Absolute Percentage Error (MAPE) between state estimation (SE) results  and power flow (PF) results.

**Arguments**
- `SE_RES::Dict`: A dictionary containing the state estimation results.
  
- `PF_RES::Dict`: A dictionary containing the power flow results.
  
- `math`: A mathematical utility or context used for processing the solutions.
  

**Returns**
- `mean_APE::Float64`: The mean absolute percentage error across all buses and terminals.
  
- `APEs::Vector{Float64}`: A vector containing the absolute percentage errors for each bus and terminal.
  

**Details**

The function compares the voltage solutions from the state estimation (`SE_RES`) and  power flow (`PF_RES`) results. It calculates the absolute percentage error (APE) for  each bus and terminal, filters out any `NaN` values, and computes the mean APE.

**Notes**
- The function assumes that the solutions in `SE_RES` and `PF_RES` are structured as  dictionaries with nested keys for buses and their respective voltages.
  
- The `dictify_solution!` function is used to process the solutions before comparison.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers._phase_letter-Tuple{Int64}' href='#Pliers._phase_letter-Tuple{Int64}'><span class="jlbinding">Pliers._phase_letter</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_phase_letter(phase::Int)
```


Convert a phase index to its corresponding letter designation.

**Arguments**
- `phase::Int`: Phase index (0=ground, 1=a, 2=b, 3=c, 4=neutral).
  

**Returns**

A single character string representing the phase.

**Throws**
- `ErrorException`: If phase is not in range 0-4.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers._pretty_diag_matrix-Tuple{Matrix}' href='#Pliers._pretty_diag_matrix-Tuple{Matrix}'><span class="jlbinding">Pliers._pretty_diag_matrix</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_pretty_diag_matrix(mat::Matrix)
```


Print a matrix as a pretty table with highlighted diagonal and off-diagonal elements.

**Arguments**
- `mat::Matrix`: The matrix to display.
  

**Description**

Displays the matrix with row indices prepended and highlights:
- Diagonal elements in cyan
  
- Off-diagonal elements in white
  
- Row labels in bold
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers._separate_phase_neutral_voltages-Tuple{Any, Any}' href='#Pliers._separate_phase_neutral_voltages-Tuple{Any, Any}'><span class="jlbinding">Pliers._separate_phase_neutral_voltages</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_separate_phase_neutral_voltages(pf_sol, bus_index)
```


Separate the phase and neutral voltages from the power flow solution for a given bus.

**Arguments**
- `pf_sol::Dict`: The power flow solution dictionary containing voltage information.
  
- `bus_index::Int`: The index of the bus for which to separate the voltages.
  

**Returns**
- `phase_voltage::Vector{ComplexF64}`: A vector containing the phase voltages (up to 3 phases).
  
- `neutral_voltage::ComplexF64`: The neutral voltage. If the neutral voltage is not present, returns 0 + 0im &quot;assuming grounded&quot;.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.add_linecode_sequences!-Tuple{Any}' href='#Pliers.add_linecode_sequences!-Tuple{Any}'><span class="jlbinding">Pliers.add_linecode_sequences!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
add_linecode_sequences!(data_eng)
```


Add linecode sequences to the `data_eng` dictionary.

**Arguments**
- `data_eng`: A dictionary containing network data.
  

This function reads linecode data from a file and adds the linecode sequences to the `data_eng` dictionary.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.add_needed_details!-Tuple{Any}' href='#Pliers.add_needed_details!-Tuple{Any}'><span class="jlbinding">Pliers.add_needed_details!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
add_needed_details!(data_eng)
```


Add additional details to the `data_eng` dictionary for lines and loads.

**Arguments**
- `data_eng`: A dictionary containing engineering data.
  

**Description**

This function adds the following details to the `data_eng` dictionary for lines:
- `"Rs (Ω)"`: The resistance of the line calculated as the product of the line length and the resistance per unit length (`rs`).
  
- `"Xs (Ω)"`: The reactance of the line calculated as the product of the line length and the reactance per unit length (`xs`).
  

For loads, the function adds the following details:
- `"PF"`: The power factor calculated as the ratio of the active power (`pd_nom`) to the apparent power (`sqrt(pd_nom^2 + qd_nom^2)`).
  
- `"S (VAh)"`: The apparent power calculated as the square root of the sum of the squares of the active power and reactive power (`sqrt(pd_nom^2 + qd_nom^2)`).
  
- `"connected_phase"`: The phase connection of the load, represented as a string (`"a"`, `"b"`, `"c"`, or `"G"`) based on the value of the `connections` field.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.add_vcomplex_vm_va!-Tuple{Any}' href='#Pliers.add_vcomplex_vm_va!-Tuple{Any}'><span class="jlbinding">Pliers.add_vcomplex_vm_va!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
add_vcomplex_vm_va!(res)
```


This function calculates the complex voltage, magnitude, and angle for each bus in the given network data.

**Arguments**
- `res`: A dictionary containing the network data.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.add_vmn_p_q-Tuple{Any, Any}' href='#Pliers.add_vmn_p_q-Tuple{Any, Any}'><span class="jlbinding">Pliers.add_vmn_p_q</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
add_vmn_p_q(math, pf_sol) -> math_meas
```


Adds voltage magnitude (VMN), active power (P), and reactive power (Q) measurements  to the given mathematical model.

**Arguments**
- `math`: The mathematical model to which the measurements will be added.
  
- `pf_sol`: The power flow solution containing the necessary data for VMN, P, and Q.
  

**Returns**
- `math_meas`: A deep copy of the input `math` with the added VMN, P, and Q measurements.
  

**Notes**

This function internally calls `_get_vmn` to add voltage magnitude measurements  and `_get_pd_qd` to add active and reactive power measurements.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.attach_Yg_ground!-Tuple{Any}' href='#Pliers.attach_Yg_ground!-Tuple{Any}'><span class="jlbinding">Pliers.attach_Yg_ground!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
attach_Yg_ground!(eng)
```


This function attaches the Yg ground to the specified transformer in the given `eng` structure.

**Arguments**
- `eng`: The structure containing the network data.
  

**Description**
- The function sets the base frequency to 50 and the default sbase to 22000.
  
- For each transformer in the `eng` structure, it finds the grounded bus.
  
- It then calculates the grnd_react value based on the line data.
  
- The function updates the `xg` and `rg` values of the grounded bus using the `grnd_react` value.
  
- Finally, it deletes the line with the `grnd_react` value from the `eng` structure.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.attach_linecodes!-Tuple{Any}' href='#Pliers.attach_linecodes!-Tuple{Any}'><span class="jlbinding">Pliers.attach_linecodes!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
attach_linecodes!(data_eng)
```


Attach linecodes to the given network data.

**Arguments**
- `data_eng`: A dictionary containing network data.
  

**Description**

This function attaches linecodes to the lines in the network data. It retrieves the linecode information from the `linecode` dictionary and assigns the corresponding values to the line dictionary.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.attach_loads!-Tuple{Any}' href='#Pliers.attach_loads!-Tuple{Any}'><span class="jlbinding">Pliers.attach_loads!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
attach_loads!(data_eng)
```


Attach loads to the specified bus in the given `data_eng` dictionary.

**Arguments**
- `data_eng`: A dictionary containing network data.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.author-Tuple{}' href='#Pliers.author-Tuple{}'><span class="jlbinding">Pliers.author</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
author()
```


Print information about the package author.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.bus_phasor!-Tuple{Makie.PolarAxis, Dict{String, Any}, Integer}' href='#Pliers.bus_phasor!-Tuple{Makie.PolarAxis, Dict{String, Any}, Integer}'><span class="jlbinding">Pliers.bus_phasor!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
bus_phasor!(ax::PolarAxis, eng::Dict{String, Any}, bus_id::String;
            linestyle=:solid, colors=[:darkred, :darkgreen, :darkblue, :black])
```


Plots voltage phasors for a given bus onto an existing `PolarAxis`.

**Arguments**
- `ax::PolarAxis`: The Makie PolarAxis to plot on.
  
- `eng::Dict{String, Any}`: The engineering data dictionary containing results.
  
- `bus_id::String`: The ID of the bus to plot phasors for.
  
- `linestyle`: The line style for the phasors (default: `:solid`).
  
- `colors`: The colors to use for the different phases (default: `[:darkred, :darkgreen, :darkblue, :black]`).
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.bus_phasor-Tuple{Dict{String, Any}, Integer}' href='#Pliers.bus_phasor-Tuple{Dict{String, Any}, Integer}'><span class="jlbinding">Pliers.bus_phasor</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
bus_phasor(eng::Dict{String, Any}, bus_id::String;
           makie_backend=WGLMakie, fig_size=(800, 800))
```


Creates a new figure and plots voltage phasors for a given bus.

**Arguments**
- `eng::Dict{String, Any}`: The engineering data dictionary containing results.
  
- `bus_id::String`: The ID of the bus to plot phasors for.
  
- `makie_backend`: The Makie backend to activate (default: `WGLMakie`).
  
- `fig_size`: The size of the figure (default: `(800, 800)`).
  

**Returns**
- `Figure`: The Makie figure containing the phasor plot.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.buses_table-Tuple{Dict{String, Any}, Function}' href='#Pliers.buses_table-Tuple{Dict{String, Any}, Function}'><span class="jlbinding">Pliers.buses_table</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
buses_table(eng::Dict{String,Any}, condition)
```


Generate and display a filtered table of buses from the given engineering data.

**Arguments**
- `eng::Dict{String,Any}`: A dictionary containing engineering data, which must include a &quot;bus&quot; key with bus information.
  
- `condition`: A function that takes a bus dictionary as input and returns a boolean indicating whether the bus meets the filtering criteria.
  

**Description**

This function extracts bus information from the provided engineering data dictionary, applies the given condition to filter the buses, and then displays the filtered buses in a formatted table. Each bus is augmented with its `bus_id` before filtering. The table includes columns for `bus_id`, `status`, `terminals`, `rg`, `xg`, and `grounded`.

**Example**

`````julia
using PowerModelsDistribution
using Pliers 
eng= Pliers.parse_file("example.dss")
buses_table(eng, bus -> bus["bus_id"] =="sourcebus")

````
or 

`````


julia buses_table(eng, bus -&gt; haskey(bus, &quot;grounded&quot;) &amp;&amp; bus[&quot;grounded&quot;]==[4]) ```


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.buses_table-Tuple{Dict{String, Any}}' href='#Pliers.buses_table-Tuple{Dict{String, Any}}'><span class="jlbinding">Pliers.buses_table</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
buses_table(eng::Dict{String, Any})
```


Generate a table summarizing the buses in the electrical network described by the dictionary `eng`.

**Arguments**
- `eng::Dict{String, Any}`: A dictionary containing various components of the electrical network.
  

**Description**

This function extracts the buses from the `eng` dictionary and creates a DataFrame with the bus ID, status, terminals, resistance to ground (rg), reactance to ground (xg), and grounding status. It then prints a formatted table of the buses.

**Example**

```julia
using PowerModelsDistribution
using Pliers 
eng= PowerModelsDistribution.parse_file("example.dss")
buses_table(eng)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.buses_table-Tuple{String}' href='#Pliers.buses_table-Tuple{String}'><span class="jlbinding">Pliers.buses_table</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Processes a DSS file and calls the function with the parsed data.

**Arguments**
- `dss::String`: The path to the DSS file to be processed.
  

**Returns**

Whatever function returns when called with parsed data.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.calc_bases_from_dict-Tuple{Dict{String, Any}}' href='#Pliers.calc_bases_from_dict-Tuple{Dict{String, Any}}'><span class="jlbinding">Pliers.calc_bases_from_dict</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



calc_bases_from_dict(data::Dict{String,Any}; return_dict=false) 

Calculate electrical base quantities from a data dictionary.

**Arguments**
- `data::Dict{String,Any}`: A dictionary containing network data. It must include:
  - `"per_unit"`: (optional) A boolean indicating whether the values are in per-unit.
    
  - `"settings"`: A nested dictionary with the following keys:
    - `"vbases_default"`: A collection where the first element has a `second` field representing the base voltage in volts.
      
    - `"voltage_scale_factor"`: A scaling factor for voltage.
      
    - `"sbase_default"`: The default base apparent power in VA.
      
    - `"power_scale_factor"`: A scaling factor for power.
      
    
  
- `return_dict::Bool=false`: If `true`, returns a dictionary of base quantities; otherwise, returns a tuple.
  

**Returns**

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
  

**Example**
- is_perunit, vbase_V, sbase_VA, Zbase_Ω, Ibase_A, vbase_ll, Ibase_A_ll, Ibase_A_ϕ = calc_bases_from_dict(data)
  
- bases = calc_bases_from_dict(data; return_dict=true)
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.calculate_vuf!-Tuple{Any, Any}' href='#Pliers.calculate_vuf!-Tuple{Any, Any}'><span class="jlbinding">Pliers.calculate_vuf!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
calculate_vuf(PF_RES, math)
```


Calculates the Voltage Unbalance Factor (VUF) for each bus in the power flow result.

VUF is defined as:     VUF = |V₋ / V₊| × 100% where:
- V₊ is the positive-sequence voltage
  
- V₋ is the negative-sequence voltage
  

It uses the phase-to-neutral complex voltages from the power flow solution, assumed to be available at each bus under keys &quot;1&quot;, &quot;2&quot;, &quot;3&quot; (phases) and &quot;4&quot; (neutral) after applying `dictify_solution!`.

Returns the modified PF_RES with VUF values added at each bus as:     PF_RES[&quot;solution&quot;][&quot;bus&quot;][bus_id][&quot;vuf&quot;] = &lt;Float64&gt;


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.convert_keys_to_symbols-Tuple{Any}' href='#Pliers.convert_keys_to_symbols-Tuple{Any}'><span class="jlbinding">Pliers.convert_keys_to_symbols</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
convert_keys_to_symbols(data)
```


Recursively convert all dictionary keys from strings to symbols.

**Arguments**
- `data`: Input data structure (Dict, Vector, or other types).
  

**Returns**
- If `data` is a Dict: A new Dict with Symbol keys and recursively converted values.
  
- If `data` is a Vector: A new Vector with recursively converted elements.
  
- Otherwise: The original value unchanged.
  

**Examples**

```julia
data = Dict("key1" => Dict("nested" => 1), "key2" => [Dict("a" => 2)])
result = convert_keys_to_symbols(data)
# result[:key1][:nested] == 1
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.create_network_graph-Tuple{Dict{String, Any}, Any}' href='#Pliers.create_network_graph-Tuple{Dict{String, Any}, Any}'><span class="jlbinding">Pliers.create_network_graph</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
create_network_graph(data::Dict{String,Any}, fallback_layout)
```


Create a MetaDiGraph representation of a power distribution network.

Automatically detects whether the input data is in ENGINEERING or MATHEMATICAL format and calls the appropriate graph creation function.

**Arguments**
- `data::Dict{String,Any}`: Network data dictionary (either ENGINEERING or MATHEMATICAL model).
  
- `fallback_layout`: Layout function to use if no coordinates are available.
  

**Returns**

The result of either `create_network_graph_eng` or `create_network_graph_math` depending on the data model type.

**See also**
- [`create_network_graph_eng`](/api#Pliers.create_network_graph_eng-Tuple{Dict{String,%20Any},%20Any})
  
- [`create_network_graph_math`](/api#Pliers.create_network_graph_math-Tuple{Dict{String,%20Any},%20Any})
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.create_network_graph_eng-Tuple{Dict{String, Any}, Any}' href='#Pliers.create_network_graph_eng-Tuple{Dict{String, Any}, Any}'><span class="jlbinding">Pliers.create_network_graph_eng</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
create_network_graph_eng(eng::Dict{String,Any}, fallback_layout) -> MetaDiGraph, Function, Dict{Symbol,Any}
```


Create a MetaDiGraph from an ENGINEERING model network.

Constructs a directed graph where vertices represent buses and edges represent lines and transformers. Bus properties are enriched with connected loads, shunts, and voltage sources.

**Arguments**
- `eng::Dict{String,Any}`: Network data in ENGINEERING format.
  
- `fallback_layout`: Layout function to use if bus coordinates are not available.
  

**Returns**

A tuple containing:
- `network_graph::MetaDiGraph`: The constructed network graph.
  
- `GraphLayout::Function`: Layout function for plotting (coordinate-based or fallback).
  
- `eng_sym::Dict{Symbol,Any}`: Copy of input data with keys converted to symbols.
  

**Examples**

```julia
network_graph, layout, eng_sym = create_network_graph_eng(eng, GraphMakie.Buchheim())
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.create_network_graph_math-Tuple{Dict{String, Any}, Any}' href='#Pliers.create_network_graph_math-Tuple{Dict{String, Any}, Any}'><span class="jlbinding">Pliers.create_network_graph_math</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
create_network_graph_math(math::Dict{String,Any}, fallback_layout) -> MetaDiGraph, Function, Dict{Symbol,Any}
```


Create a MetaDiGraph from a MATHEMATICAL model network.

Constructs a directed graph where vertices represent buses and edges represent branches. Bus properties are enriched with connected loads, shunts, and generators.

**Arguments**
- `math::Dict{String,Any}`: Network data in MATHEMATICAL format.
  
- `fallback_layout`: Layout function to use if bus coordinates are not available.
  

**Returns**

A tuple containing:
- `network_graph::MetaDiGraph`: The constructed network graph.
  
- `GraphLayout::Function`: Layout function for plotting (coordinate-based or fallback).
  
- `math_sym::Dict{Symbol,Any}`: Copy of input data with keys converted to symbols.
  

**Examples**

```julia
network_graph, layout, math_sym = create_network_graph_math(math, GraphMakie.Buchheim())
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.default_index_value-Tuple{Integer, Symbol, Set}' href='#Pliers.default_index_value-Tuple{Integer, Symbol, Set}'><span class="jlbinding">Pliers.default_index_value</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
default_index_value(v, prop, index_values; exclude=nothing)
```


Provides a default index value for a vertex if no value currently exists. The default is a string: &quot;$prop$i&quot; where `prop` is the property name and `i` is the vertex number. If some other vertex already has this name, a randomized string is generated (though the way it is generated is deterministic).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.dictify_solution!-Tuple{Dict{String, Any}, Dict{String, Any}}' href='#Pliers.dictify_solution!-Tuple{Dict{String, Any}, Dict{String, Any}}'><span class="jlbinding">Pliers.dictify_solution!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
dictify_solution!(pf_sol::Dict{String, Any}, math::Dict{String, Any}; formulation = "IVR")
```


Transforms and organizes the solution data from a power flow computation into a structured dictionary format.

**Arguments**
- `pf_sol::Dict{String, Any}`: A dictionary containing the power flow solution data.
  
- `math::Dict{String, Any}`: A dictionary containing the mathematical model data.
  
- `formulation::String` (optional): Specifies the formulation type to be used. Defaults to `"IVR"`.
  

**Description**

This function modifies the `pf_sol` dictionary in-place by calling helper functions to process and structure data for buses, loads, branches, and generators. Each helper function is responsible for handling a specific component of the power flow solution.

**Notes**
- The function assumes that the helper functions `solution_dictify_buses!`, `solution_dictify_loads!`, `solution_dictify_branches!`, and `solution_dictify_gens!` are defined and properly handle their respective components.
  
- The `formulation` parameter allows customization of the solution processing based on the formulation type.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.diff_vectors-Tuple{Vector{Float64}, Vector{Float64}}' href='#Pliers.diff_vectors-Tuple{Vector{Float64}, Vector{Float64}}'><span class="jlbinding">Pliers.diff_vectors</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
diff_vectors(vec1::Vector{Float64}, vec2::Vector{Float64})
```


Prints the difference between two vectors element-wise.

**Arguments**
- `vec1::Vector{Float64}`: The first vector.
  
- `vec2::Vector{Float64}`: The second vector.
  

**Example**

diff_vectors([1.0, 2.0, 3.0], [1.0, 2.0, 4.0])


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.eng_report-Tuple{Dict{String, Any}}' href='#Pliers.eng_report-Tuple{Dict{String, Any}}'><span class="jlbinding">Pliers.eng_report</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
eng_report(eng::Dict{String, Any})
```


Generate a report for the electrical network described by the dictionary `eng`.

**Arguments**
- `eng::Dict{String, Any}`: A dictionary containing various components of the electrical network.
  

**Description**

This function extracts various components from the `eng` dictionary, such as buses, lines, linecodes, loads, voltage sources, time series data, conductor IDs, name, settings, files, and data model. It then prints a formatted report summarizing the contents of the network, including the number of each component present.

**Example**

```julia
using PowerModelsDistribution
using Pliers 
eng= PowerModelsDistribution.parse_file("example.dss")
eng_report(eng)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.error_text-Tuple{String}' href='#Pliers.error_text-Tuple{String}'><span class="jlbinding">Pliers.error_text</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
error_text(message::String)
```


Print an error message with italic red formatting.

**Arguments**
- `message::String`: The error message to display.
  

**Examples**

```julia
error_text("This is an error")
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.extra_keys-Tuple{Dict{String, Any}, Any}' href='#Pliers.extra_keys-Tuple{Dict{String, Any}, Any}'><span class="jlbinding">Pliers.extra_keys</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
extra_keys(eng_data::Dict{String, Any}, expected_keys)
```


Checks if there are extra keys in the KeySet of the dictionary `eng_data` that are not in the list `keys` and prints them.

**Arguments**
- `eng_data::Dict{String, Any}`: A dictionary containing the data to be checked.
  
- `expected_keys`: A list of keys that are expected to be in the dictionary.
  
- `show_keys=false`: A boolean indicating whether to print the expected keys, the keys in the dictionary, and the extra keys.
  

**Description**

This function compares the keys in the dictionary `eng_data` with the list of `expected_keys` and prints a warning message if there are extra keys in the dictionary that are not in the list. If `show_keys` is set to `true`, the expected keys, the keys in the dictionary, and the extra keys are printed.

**Example**

```julia
using PowerModelsDistribution
using Pliers
eng= PowerModelsDistribution.parse_file("example.dss")
extra_keys(eng, ["bus", "line", "linecode", "load", "voltage_source", "time_series", "conductor_ids", "name", "settings", "files", "data_model"])
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.fluff_bus_voltages!-Tuple{Dict{String, Any}}' href='#Pliers.fluff_bus_voltages!-Tuple{Dict{String, Any}}'><span class="jlbinding">Pliers.fluff_bus_voltages!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
fluff_bus_voltages(PF_Res::Dict{String, Any})
```


This function is used to get all possible voltage forms from the output of the OPF voltage variables, so it processes the bus voltages in the given Power Flow Results dictionary `PF_Res`. For each bus in `PF_Res["bus"]`, it checks if the bus has real (`vr`) and imaginary (`vi`) voltage components. If both components are present, it calculates the complex voltage `V`, the voltage magnitude `vm`, and the voltage angle `va` for the bus. A warning is issued indicating that only `vr` and `vi` can be processed.

**Arguments**
- `PF_Res::Dict{String, Any}`: A dictionary containing the Power Flow Results solution, including bus voltage information.
  

!!! You need to pass the [&quot;solution&quot;] key of the Power Flow Results dictionary.

**Modifies**
- Adds the following keys to each bus dictionary if `vr` and `vi` are present:
  - `V`: The complex voltage calculated as `vr + vi*im`.
    
  - `vm`: The magnitude of the complex voltage.
    
  - `va`: The angle of the complex voltage.
    
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.get_graph_edge-Tuple{Any, Any}' href='#Pliers.get_graph_edge-Tuple{Any, Any}'><span class="jlbinding">Pliers.get_graph_edge</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_graph_edge(G, edge_id)
get_graph_edge(G, edge_id, key)
```


Get properties of an edge (line/branch) from the network graph.

**Arguments**
- `G::MetaDiGraph`: The network graph.
  
- `edge_id`: Edge identifier (line ID for ENGINEERING, branch ID for MATHEMATICAL).
  
- `key`: (Optional) Specific property key to retrieve.
  

**Returns**
- Without `key`: Dictionary of all edge properties.
  
- With `key`: Value of the specified property.
  

**Examples**

```julia
props = get_graph_edge(G, "line1")
impedance = get_graph_edge(G, "line1", "rs")
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.get_graph_node-Tuple{Any, Any}' href='#Pliers.get_graph_node-Tuple{Any, Any}'><span class="jlbinding">Pliers.get_graph_node</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_graph_node(G, node)
get_graph_node(G, node, key)
```


Get properties of a node (bus) from the network graph.

**Arguments**
- `G::MetaDiGraph`: The network graph.
  
- `node`: Node identifier (bus ID).
  
- `key`: (Optional) Specific property key to retrieve.
  

**Returns**
- Without `key`: Dictionary of all node properties.
  
- With `key`: Value of the specified property.
  

**Examples**

```julia
props = get_graph_node(G, "bus1")
voltage = get_graph_node(G, "bus1", "vm")
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.get_sequence_components-Tuple{Matrix}' href='#Pliers.get_sequence_components-Tuple{Matrix}'><span class="jlbinding">Pliers.get_sequence_components</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
get_sequence_components(m_phase::Matrix{ComplexF64})
```


Calculates the zero, positive, and negative sequence components of a 3x3 complex matrix representing phase quantities (e.g., impedance or admittance).

**Arguments**
- `m_phase::Matrix{ComplexF64}`: The input 3x3 complex matrix in phase coordinates.
  

**Returns**
- `Tuple{Matrix{ComplexF64}, ComplexF64, ComplexF64, ComplexF64}`: A tuple containing:
  1. The full 3x3 sequence matrix (M_seq).
    
  2. The zero sequence component (diagonal element M_seq[1, 1]).
    
  3. The positive sequence component (diagonal element M_seq[2, 2]).
    
  4. The negative sequence component (diagonal element M_seq[3, 3]).
    
  

**Throws**
- `DimensionMismatch`: If the input matrix is not 3x3.
  

**Usage**

Assume &#39;M_red&#39; is the 3x3 matrix obtained from the previous Kron reduction step Let&#39;s create a sample symmetrical 3x3 matrix for demonstration: Z_phase = ComplexF64[     10+5im   2+1im   2+1im;     2+1im   10+5im  2+1im;     2+1im   2+1im   10+5im ] This represents a balanced system where off-diagonals are equal.

Or use a more general (unbalanced) example: Z_phase = ComplexF64[      10+5im   1+0.5im  2+1im;      1+0.5im  12+6im   3+1.5im;      2+1im    3+1.5im  11+5.5im ]

try     # Calculate sequence components     Z_seq_matrix, Z0, Z1, Z2 = calculate_sequence_components(Z_phase)

```julia
println("Original Phase Matrix Z_phase:")
display(Z_phase)

println("
```


Sequence Matrix Z_seq:&quot;)     display(round.(Z_seq_matrix; digits=4)) # Round for cleaner display

```julia
println("
```


Sequence Components:&quot;)     println(&quot;Zero Sequence (Z0): &quot;, round(Z0; digits=4))     println(&quot;Positive Sequence (Z1): &quot;, round(Z1; digits=4))     println(&quot;Negative Sequence (Z2): &quot;, round(Z2; digits=4))

catch e     println(&quot;Error: &quot;, e) end


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.header-Tuple{String}' href='#Pliers.header-Tuple{String}'><span class="jlbinding">Pliers.header</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
header(text::String)
```


Print a header text with bold, underlined blue formatting.

**Arguments**
- `text::String`: The text to display as a header.
  

**Examples**

```julia
header("Section Title")
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.isolate_feeder-Tuple{Any, Any}' href='#Pliers.isolate_feeder-Tuple{Any, Any}'><span class="jlbinding">Pliers.isolate_feeder</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



isolate_feeder(eng, feederno)

Given a network data structure `eng` and a feeder number `feederno`, this function creates a deep copy of the `eng` structure and isolates the specified feeder.

**Arguments**
- `eng`: A network data structure.
  
- `feederno`: The feeder number to isolate.
  

**Returns**
- `eng_feeder`: A deep copy of the `eng` structure with the specified feeder isolated.
  

**Description**

The function isolates the specified feeder by performing the following steps:
1. Creates a deep copy of the `eng` structure.
  
2. Sets the voltage source bus to the TrLVBus of the specified feeder.
  
3. Retrieves the transformer ID(s) associated with the TrLVBus.
  
4. Retrieves the feeder head bus and the circuit breaker buses of the specified feeder.
  
5. Constructs a list of member lines by traversing the network starting from the feeder head bus.
  
6. Filters out any empty member lines and flattens the list.
  
7. Appends the transformer ID(s), feeder number, and circuit breaker number to the list of member lines.
  
8. Removes irrelevant data from the `eng_feeder` structure based on the transformer ID(s), feeder number, and member lines.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.linecodes_table-Tuple{String}' href='#Pliers.linecodes_table-Tuple{String}'><span class="jlbinding">Pliers.linecodes_table</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Processes a DSS file and calls the function with the parsed data.

**Arguments**
- `dss::String`: The path to the DSS file to be processed.
  

**Returns**

Whatever function returns when called with parsed data.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.lines_table-Tuple{Dict{String, Any}, Function}' href='#Pliers.lines_table-Tuple{Dict{String, Any}, Function}'><span class="jlbinding">Pliers.lines_table</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lines_table(eng::Dict{String,Any}, condition)
```


Generate and display a filtered table of lines from the given engineering data.

**Arguments**
- `eng::Dict{String,Any}`: A dictionary containing engineering data, which must include a &quot;line&quot; key with line information.
  
- `condition`: A function that takes a line dictionary as input and returns a boolean indicating whether the line meets the filtering criteria.
  

**Description**

This function extracts line information from the provided engineering data dictionary, applies the given condition to filter the lines, and then displays the filtered lines in a formatted table. Each line is augmented with its `line_id` before filtering. The table includes columns for `line_id`, `source_id`, `f_bus`, `f_connections`, `t_bus`, `t_connections`, `length`, `linecode`, and `status`.

**Example**

```julia
using PowerModelsDistribution
using Pliers
eng= PowerModelsDistribution.parse_file("example.dss")
lines_table(eng, line -> line["length"] > 0.75)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.lines_table-Tuple{Dict{String, Any}}' href='#Pliers.lines_table-Tuple{Dict{String, Any}}'><span class="jlbinding">Pliers.lines_table</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lines_table(eng::Dict{String, Any})
```


Generate a table summarizing the lines in the electrical network described by the dictionary `eng`.

**Arguments**
- `eng::Dict{String, Any}`: A dictionary containing various components of the electrical network.
  

**Description**

This function extracts the lines from the `eng` dictionary and creates a DataFrame with the line ID, status, from bus, to bus, length, resistance, reactance, and linecode. It then prints a formatted table of the lines.

**Example**

```julia
using PowerModelsDistribution
using Pliers
eng= PowerModelsDistribution.parse_file("example.dss")
lines_table(eng)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.lines_table-Tuple{String}' href='#Pliers.lines_table-Tuple{String}'><span class="jlbinding">Pliers.lines_table</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Processes a DSS file and calls the function with the parsed data.

**Arguments**
- `dss::String`: The path to the DSS file to be processed.
  

**Returns**

Whatever function returns when called with parsed data.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.list_files-Tuple{Any}' href='#Pliers.list_files-Tuple{Any}'><span class="jlbinding">Pliers.list_files</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
list_files(directory)
```


List all files recursively in a directory.

**Arguments**
- `directory`: Root directory to list files from.
  

**Returns**

A vector of full paths to all files in the directory tree.

**Examples**

```julia
all_files = list_files("/path/to/dir")
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.loads_table-Tuple{Dict{String, Any}, Function}' href='#Pliers.loads_table-Tuple{Dict{String, Any}, Function}'><span class="jlbinding">Pliers.loads_table</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
loads_table(eng::Dict{String,Any}, condition)
```


Generate and display a filtered table of loads from the given engineering data.

**Arguments**
- `eng::Dict{String,Any}`: A dictionary containing engineering data, which must include a &quot;load&quot; key with load information.
  
- `condition`: A function that takes a load dictionary as input and returns a boolean indicating whether the load meets the filtering criteria.
  

**Description**

This function extracts load information from the provided engineering data dictionary, applies the given condition to filter the loads, and then displays the filtered loads in a formatted table. Each load is augmented with its `load_id` before filtering. The table includes columns for `load_id`, `source_id`, `bus`, `connections`, `vm_nom`, `pd_nom`, `qd_nom`, `configuration`, `model`, `dispatchable`, and `status`.

**Example**
- filter the loads by the value of `pd_nom`:   `julia   loads_table(eng, load -> load["pd_nom"] > [0.33])`
  
- filter the loads by the phase connectivity   `julia   loads_table(eng, load -> load["connections"] == [1, 4])`
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.loads_table-Tuple{Dict{String, Any}}' href='#Pliers.loads_table-Tuple{Dict{String, Any}}'><span class="jlbinding">Pliers.loads_table</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
loads_table(eng::Dict{String, Any})
```


Generate a table summarizing the loads in the electrical network described by the dictionary `eng`.

**Arguments**
- `eng::Dict{String, Any}`: A dictionary containing various components of the electrical network.
  

**Description**

This function extracts the loads from the `eng` dictionary and creates a DataFrame with the load ID, status, bus, phases, kw, kvar, and kva. It then prints a formatted table of the loads.

**Example**

```julia
using PowerModelsDistribution
using Pliers
eng= PowerModelsDistribution.parse_file("example.dss")
loads_table(eng)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.loads_table-Tuple{String}' href='#Pliers.loads_table-Tuple{String}'><span class="jlbinding">Pliers.loads_table</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Processes a DSS file and calls the function with the parsed data.

**Arguments**
- `dss::String`: The path to the DSS file to be processed.
  

**Returns**

Whatever function returns when called with parsed data.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.math_transformers_table-Tuple{Dict{String, Any}, Function}' href='#Pliers.math_transformers_table-Tuple{Dict{String, Any}, Function}'><span class="jlbinding">Pliers.math_transformers_table</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
math_transformers_table(math::Dict{String, Any}, condition::Function)
```


Generate and display a filtered table of transformers from the given mathematical model data.

**Arguments**
- `math::Dict{String,Any}`: A dictionary containing mathematical model data, which must include a &quot;transformer&quot; key with transformer information.
  
- `condition`: A function that takes a transformer dictionary as input and returns a boolean indicating whether the transformer meets the filtering criteria.
  

**Description**

This function extracts transformer information from the provided mathematical model dictionary, applies the given condition to filter the transformers, and then displays the filtered transformers in a formatted table.

**Example**

```julia
using PowerModelsDistribution
using Pliers
eng = PowerModelsDistribution.parse_file("example.dss")
math = PowerModelsDistribution.transform_data_model(eng)
math_transformers_table(math, transformer -> transformer["f_bus"] == 1)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.math_transformers_table-Tuple{Dict{String, Any}}' href='#Pliers.math_transformers_table-Tuple{Dict{String, Any}}'><span class="jlbinding">Pliers.math_transformers_table</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
math_transformers_table(math::Dict{String, Any})
```


Generate a table summarizing the transformers in the mathematical model described by the dictionary `math`.

**Arguments**
- `math::Dict{String, Any}`: A dictionary containing various components of the mathematical model.
  

**Description**

This function extracts the transformers from the `math` dictionary and creates a DataFrame with the transformer index, name, source ID, from/to buses, connections, configuration, tap settings, voltage bases, polarity, power/current limits, and status. It then prints a formatted table of the transformers.

**Example**

```julia
using PowerModelsDistribution
using Pliers 
eng = PowerModelsDistribution.parse_file("example.dss")
math = PowerModelsDistribution.transform_data_model(eng)
math_transformers_table(math)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.merge_bus_results!-Tuple{Any, Any}' href='#Pliers.merge_bus_results!-Tuple{Any, Any}'><span class="jlbinding">Pliers.merge_bus_results!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
merge_bus_results!(res_mn, mn_Feeder_math)
```


Merge the bus results from `res_mn` into `mn_Feeder_math`.

**Arguments**
- `res_mn`: A dictionary containing the bus results.
  
- `mn_Feeder_math`: A dictionary representing the network data.
  

**Description**

This function merges the bus results from `res_mn` into the `mn_Feeder_math` dictionary. It updates the `"vi"`, `"va"`, and `"vcomplex"` values for each bus in the network.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.network_graph_map_plot-Tuple{MetaGraphs.MetaDiGraph, Function}' href='#Pliers.network_graph_map_plot-Tuple{MetaGraphs.MetaDiGraph, Function}'><span class="jlbinding">Pliers.network_graph_map_plot</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
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
```


Plots a network graph on a map using the specified layout and visual properties.

**Arguments**
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
  

**Returns**
- A map object with the plotted network graph.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.network_graph_plot-Tuple{MetaGraphs.MetaDiGraph}' href='#Pliers.network_graph_plot-Tuple{MetaGraphs.MetaDiGraph}'><span class="jlbinding">Pliers.network_graph_plot</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
network_graph_plot(network_graph::MetaDiGraph; layout=GraphMakie.Buchheim(), figure_size=(1000, 1200), node_color=:black, node_size=automatic, node_marker=automatic, node_strokewidth=automatic, show_node_labels=false, show_edge_labels=false, edge_color=:black, elabels_color=:black, elabels_fontsize=10, tangents=((0,-1),(0,-1)), arrow_show=false, arrow_marker='➤', arrow_size=50, arrow_shift=0.5)
```


Plots the given network graph using the specified layout and styling options.

**Arguments**
- `network_graph::MetaDiGraph`: The network graph to be plotted.
  

**Keyword Arguments**
- `layout`: The layout algorithm to use for positioning the nodes (default: `GraphMakie.Buchheim()`).
  
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
  

**Returns**
- A plot object representing the network graph.
  

**Example**

```julia
f, ax, p = network_graph_plot(network_graph; layout=GraphMakie.Buchheim(), figure_size=(1000, 1200), node_color=:black, node_size=automatic, node_marker=automatic, node_strokewidth=automatic, show_node_labels=false, show_edge_labels=false, edge_color=:black, elabels_color=:black, elabels_fontsize=10, tangents=((0,-1),(0,-1)), arrow_show=false, arrow_marker='➤', arrow_size=50, arrow_shift=0.5)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.plot_network-Tuple{Any}' href='#Pliers.plot_network-Tuple{Any}'><span class="jlbinding">Pliers.plot_network</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_network(eng; plot_bus=true, overlayNetwork=false, save_fig=false, plotname="plot.png"::String, save_json=false, json_filename="Networkfile.json"::String)
```


This function plots the network based on the given network data.

**Arguments**
- `eng`: A dictionary containing the network data.
  
- `plot_bus`: A boolean indicating whether to plot the buses. Default is `true`.
  
- `overlayNetwork`: A boolean indicating whether to overlay the network on an existing plot. Default is `false`.
  
- `save_fig`: A boolean indicating whether to save the plot as an image file. Default is `false`.
  
- `plotname`: A string specifying the name of the image file to be saved. Default is &quot;plot.png&quot;.
  
- `save_json`: A boolean indicating whether to save the network data as a JSON file. Default is `false`.
  
- `json_filename`: A string specifying the name of the JSON file to be saved. Default is &quot;Networkfile.json&quot;.
  

**Returns**
- `p`: The plot object representing the network plot.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.plot_network_coords-Tuple{Dict{String, Any}}' href='#Pliers.plot_network_coords-Tuple{Dict{String, Any}}'><span class="jlbinding">Pliers.plot_network_coords</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_network_coords(eng::Dict{String,Any}; show_node_labels=false, show_edge_labels=false, fallback_layout=GraphMakie.Buchheim(), kwargs...)
```


Plot a network graph with optional node and edge labels.

**Arguments**
- `eng::Dict{String,Any}`: The engineering data used to create the network graph.
  
- `show_node_labels::Bool`: Whether to display labels for the nodes. Default is `false`.
  
- `show_edge_labels::Bool`: Whether to display labels for the edges. Default is `false`.
  
- `fallback_layout`: The layout algorithm to use if no specific layout is provided. Default is `GraphMakie.Buchheim()`.
  
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
  

**Returns**
- A plot of the network graph with the specified layout and labels.
  

**Example**

```julia
plot_network_coords(eng)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.plot_network_tree-Tuple{Dict{String, Any}}' href='#Pliers.plot_network_tree-Tuple{Dict{String, Any}}'><span class="jlbinding">Pliers.plot_network_tree</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_network_tree(eng::Dict{String,Any}; makie_backend=WGLMakie)
```


Plots a network tree based on the given engineering model dictionary `eng`.

**Arguments**
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
  

**Returns**
- A tuple `(f, ax, p)` where `f` is the figure, `ax` is the axis, and `p` is the plot object.
  

**Details**
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
  

**Errors**
- Throws an error if `sourcebus` is not found in the bus data.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.plot_network_tree-Tuple{String}' href='#Pliers.plot_network_tree-Tuple{String}'><span class="jlbinding">Pliers.plot_network_tree</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_network_tree(dss::String; kwargs...)
```


Plots a network tree from a given DSS file.

**Arguments**
- `dss::String`: The path to the DSS file containing the network data.
  
- `makie_backend`: The Makie backend to use for plotting. Defaults to `WGLMakie`.
  

**Description**

This function parses the DSS file to create an engineering model of the network, transforms loops in the model, and converts keys to symbols. It then constructs a network graph using `MetaDiGraph`, adds vertices for each bus, and sets the `sourcebus` as the root. Edges are added based on the `f_bus` and `t_bus` indices of the lines. Finally, it plots the network graph using `graphplot` with a specified layout and labels for nodes and edges.

**Returns**

A plot of the network graph.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.polarize-Tuple{Complex}' href='#Pliers.polarize-Tuple{Complex}'><span class="jlbinding">Pliers.polarize</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
polarize(complex::Complex; scale=1.0, rounding_digits=4)
```


Convert a complex number to polar notation string format.

**Arguments**
- `complex::Complex`: The complex number to convert.
  

**Keyword Arguments**
- `scale::Float64`: Scaling factor for the magnitude (default: 1.0).
  
- `rounding_digits::Int`: Number of decimal places for rounding (default: 4).
  

**Returns**

A string in the format &quot;magnitude ∠ angle_degrees&quot;.

**Examples**

```julia
polarize(1 + 1im)  # "1.4142 ∠ 45.0"
polarize(1 + 1im; scale=1000)  # "1414.2136 ∠ 45.0"
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.read_data-Tuple{String}' href='#Pliers.read_data-Tuple{String}'><span class="jlbinding">Pliers.read_data</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
read_data(path::String)
```


Reads and deserializes data from the specified `path`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.remove_virtual_bus!-Tuple{Any}' href='#Pliers.remove_virtual_bus!-Tuple{Any}'><span class="jlbinding">Pliers.remove_virtual_bus!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
remove_virtual_bus!(math_b)
```


Remove virtual bus, branch, and generator from the given `math` dictionary.

**Arguments**
- `math_b::Dict`: A dictionary containing information about buses, branches, and generators.
  

**Example**


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.rm_spanish_transformer!-Tuple{Any}' href='#Pliers.rm_spanish_transformer!-Tuple{Any}'><span class="jlbinding">Pliers.rm_spanish_transformer!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
rm_spanish_transformer!(data_eng)
```


Remove the Spanish transformer from the given `data_eng` dictionary.

**Arguments**
- `data_eng`: A dictionary containing the network data.
  

**Description**

This function removes the Spanish transformer from the `data_eng` dictionary. It updates the voltage source parameters, deletes the transformer, deletes the lines connected before the transformer, updates the bus settings, and resets the terminals and grounded status of the buses.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.save_data-Tuple{Any, String}' href='#Pliers.save_data-Tuple{Any, String}'><span class="jlbinding">Pliers.save_data</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
save_data(data, path::String)
```


Serializes and saves `data` to the specified `path`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.search_directories-Tuple{Any, Any}' href='#Pliers.search_directories-Tuple{Any, Any}'><span class="jlbinding">Pliers.search_directories</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
search_directories(directory, directory_name)
```


Search recursively for directories containing a specific string in their name.

**Arguments**
- `directory`: Root directory to start the search.
  
- `directory_name`: String to search for in directory names.
  

**Returns**

A vector of full paths to matching directories.

**Examples**

```julia
dirs = search_directories("/path/to/dir", "test")
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.search_files-Tuple{Any, Any}' href='#Pliers.search_files-Tuple{Any, Any}'><span class="jlbinding">Pliers.search_files</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
search_files(directory, file_name)
```


Search recursively for files containing a specific string in their name.

**Arguments**
- `directory`: Root directory to start the search.
  
- `file_name`: String to search for in file names.
  

**Returns**

A vector of full paths to matching files.

**Examples**

```julia
files = search_files("/path/to/dir", "test")
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.set_indexing_prop!-Tuple{MetaGraphs.AbstractMetaGraph, Symbol}' href='#Pliers.set_indexing_prop!-Tuple{MetaGraphs.AbstractMetaGraph, Symbol}'><span class="jlbinding">Pliers.set_indexing_prop!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
set_indexing_prop!(g, prop)
set_indexing_prop!(g, v, prop, val)
```


Make property `prop` into an indexing property. If any values for this property are already set, each vertex must have unique values. Optionally, set the index `val` for vertex `v`. Any vertices without values will be set to a default (&quot;(prop)(v)&quot;).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.set_journal_theme-Tuple{}' href='#Pliers.set_journal_theme-Tuple{}'><span class="jlbinding">Pliers.set_journal_theme</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
set_journal_theme(; fontsize=nothing)
```


Set a Makie theme suitable for journal publications (IEEE style).

This function configures a Makie theme with appropriate fonts, colors, and styling for publication-quality figures. The theme uses TeX Gyre Termes fonts when available, falling back to Computer Modern if not found.

**Keyword Arguments**
- `fontsize`: Optional font size override. Defaults to 8pt if not specified.
  

**Returns**

A tuple containing unit conversion factors:
- `inch::Float64`: Pixels per inch (96)
  
- `pt::Float64`: Pixels per point (4/3)
  
- `cm::Float64`: Pixels per centimeter
  
- `ieeecolumn::Float64`: Width of a single IEEE column in pixels (3.5 inches)
  
- `ieee2column::Float64`: Width of a double IEEE column in pixels (7.16 inches)
  

**Examples**

```julia
inch, pt, cm, ieeecolumn, ieee2column = set_journal_theme()
fig = Figure(size=(ieeecolumn, ieeecolumn))
```


**Notes**

IEEE Figure Sizes:
- One column width: 3.5 inches, 88.9 mm, or 21 picas
  
- Two columns width: 7.16 inches, 182 mm, or 43 picas
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.show_example-Tuple{Dict}' href='#Pliers.show_example-Tuple{Dict}'><span class="jlbinding">Pliers.show_example</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
show_example(dict::Dict)
```


Return the value of the first entry in a dictionary.

**Arguments**
- `dict::Dict`: A dictionary to inspect.
  

**Returns**

The value associated with the first key in the dictionary.

**Examples**

```julia
d = Dict("a" => 1, "b" => 2)
show_example(d)  # Returns 1 or 2 depending on iteration order
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.show_transformer_math_components-Tuple{Any}' href='#Pliers.show_transformer_math_components-Tuple{Any}'><span class="jlbinding">Pliers.show_transformer_math_components</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
show_transformer_math_components(math; suppress_print::Bool=false)
```


Extract and display transformer mathematical components from a power system model.

**Arguments**
- `math::Dict`: A dictionary containing the mathematical model representation with  &quot;map&quot; and element type keys (e.g., &quot;bus&quot;, &quot;branch&quot;, &quot;transformer&quot;)
  
- `suppress_print::Bool=false`: If `true`, suppresses console output and only returns results
  

**Returns**
- `Dict{String, Any}`: A nested dictionary organized by transformer name, containing:
  - `"buses"`: Dictionary of bus elements mapped to this transformer
    
  - `"branches"`: Dictionary of branch elements mapped to this transformer
    
  - `"transformers"`: Dictionary of transformer elements mapped to this transformer
    
  

**Description**

This function identifies all transformers in the mathematical model that were converted  from engineering model representation and extracts the associated mathematical components  (buses, branches, transformers). Each transformer is identified by the presence of an  &quot;unmap_function&quot; equal to &quot;_map_math2eng_transformer!&quot; in the mapping structure.

If `suppress_print=false`, displays formatted console output showing the transformation  process and component details.

**Example**

mapped_transformers = show_transformer_math_components(math, suppress_print=false)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.sub_header-Tuple{String}' href='#Pliers.sub_header-Tuple{String}'><span class="jlbinding">Pliers.sub_header</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sub_header(text::String)
```


Print a sub-header text with bold, italic blue formatting.

**Arguments**
- `text::String`: The text to display as a sub-header.
  

**Examples**

```julia
sub_header("Sub-section Title")
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.sub_sub_header-Tuple{String}' href='#Pliers.sub_sub_header-Tuple{String}'><span class="jlbinding">Pliers.sub_sub_header</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sub_sub_header(text::String)
```


Print a sub-sub-header text with bold magenta formatting.

**Arguments**
- `text::String`: The text to display as a sub-sub-header.
  

**Examples**

```julia
sub_sub_header("Minor Section Title")
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.transformers_table-Tuple{Dict{String, Any}, Function}' href='#Pliers.transformers_table-Tuple{Dict{String, Any}, Function}'><span class="jlbinding">Pliers.transformers_table</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
transformers_table(eng::Dict{String,Any}, condition)
```


Generate and display a filtered table of transformers from the given engineering data.

**Arguments**
- `eng::Dict{String,Any}`: A dictionary containing engineering data, which must include a &quot;transformer&quot; key with transformer information.
  
- `condition`: A function that takes a transformer dictionary as input and returns a boolean indicating whether the transformer meets the filtering criteria.
  

**Description**

This function extracts transformer information from the provided engineering data dictionary, applies the given condition to filter the transformers, and then displays the filtered transformers in a formatted table. Each transformer is augmented with its `transformer_id` before filtering.

**Example**

```julia
using PowerModelsDistribution
using Pliers
eng= PowerModelsDistribution.parse_file("example.dss")
transformers_table(eng, transformer -> transformer["vm_nom"][1] > 10.0)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.transformers_table-Tuple{Dict{String, Any}}' href='#Pliers.transformers_table-Tuple{Dict{String, Any}}'><span class="jlbinding">Pliers.transformers_table</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
transformers_table(eng::Dict{String, Any})
```


Generate a table summarizing the transformers in the electrical network described by the dictionary `eng`.

**Arguments**
- `eng::Dict{String, Any}`: A dictionary containing various components of the electrical network.
  

**Description**

This function extracts the transformers from the `eng` dictionary and creates a DataFrame with the transformer ID, name, source ID, buses, connections, nominal voltage (vm_nom), nominal power (sm_nom), configuration, short-circuit reactance (xsc), winding resistance (rw), no-load loss, magnetizing current (cmag), tap settings (tm_set), fixed taps (tm_fix), polarity, and status. It then prints a formatted table of the transformers.

**Example**

```julia
using PowerModelsDistribution
using Pliers 
eng= PowerModelsDistribution.parse_file("example.dss")
transformers_table(eng)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.transformers_table-Tuple{String}' href='#Pliers.transformers_table-Tuple{String}'><span class="jlbinding">Pliers.transformers_table</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Processes a DSS file and calls the function with the parsed data.

**Arguments**
- `dss::String`: The path to the DSS file to be processed.
  

**Returns**

Whatever function returns when called with parsed data.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Pliers.warning_text-Tuple{String}' href='#Pliers.warning_text-Tuple{String}'><span class="jlbinding">Pliers.warning_text</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
warning_text(message::String)
```


Print a warning message with italic yellow formatting.

**Arguments**
- `message::String`: The warning message to display.
  

**Examples**

```julia
warning_text("This is a warning")
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## PMDUtils Sub-module {#PMDUtils-Sub-module}

Re-exports PMD-related utility functions. These functions are also available directly from the main Pliers module.
<details class='jldocstring custom-block' open>
<summary><a id='Pliers.PMDUtils' href='#Pliers.PMDUtils'><span class="jlbinding">Pliers.PMDUtils</span></a> <Badge type="info" class="jlObjectType jlModule" text="Module" /></summary>



```julia
PMDUtils
```


Internal sub-module providing utility functions for PowerModelsDistribution (PMD) workflows.

This module re-exports functions from the main Pliers module for:
- Processing power flow solutions (voltage fluffing, dictification)
  
- Impedance calculations (Kron reduction, sequence components)
  
- Network data manipulation
  
- Result analysis and transformation
  
- Engineering and mathematical model exploration
  

See the main Pliers module for function documentation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## PMDSEUtils Sub-module {#PMDSEUtils-Sub-module}

Re-exports PMDSE-related utility functions. These functions are also available directly from the main Pliers module.
<details class='jldocstring custom-block' open>
<summary><a id='Pliers.PMDSEUtils' href='#Pliers.PMDSEUtils'><span class="jlbinding">Pliers.PMDSEUtils</span></a> <Badge type="info" class="jlObjectType jlModule" text="Module" /></summary>



```julia
PMDSEUtils
```


Internal sub-module providing utility functions for PowerModelsDistributionStateEstimation (PMDSE) workflows.

This module re-exports functions from the main Pliers module for:
- State estimation result visualization
  
- Measurement residual analysis
  
- Measurement data processing and writing
  

See the main Pliers module for function documentation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## PMDPlotting Sub-module {#PMDPlotting-Sub-module}

Re-exports network visualization and plotting functions. These functions are also available directly from the main Pliers module.
<details class='jldocstring custom-block' open>
<summary><a id='Pliers.PMDPlotting' href='#Pliers.PMDPlotting'><span class="jlbinding">Pliers.PMDPlotting</span></a> <Badge type="info" class="jlObjectType jlModule" text="Module" /></summary>



```julia
PMDPlotting
```


Internal sub-module providing plotting functions for power distribution network visualization.

This module re-exports functions from the main Pliers module for:
- Network tree visualization
  
- Coordinate-based network plotting
  
- Geographic map overlays
  
- Bus phasor diagram plotting
  

See the main Pliers module for function documentation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/MohamedNumair/Pliers.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

