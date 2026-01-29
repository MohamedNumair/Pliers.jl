"""
    PMDUtils

Internal sub-module providing utility functions for PowerModelsDistribution (PMD) workflows.

This module includes functions for:
- Processing power flow solutions (voltage fluffing, dictification)
- Impedance calculations (Kron reduction, sequence components)
- Network data manipulation
- Result analysis and transformation
- Engineering and mathematical model exploration

# Exported Functions

## Solution Processing
- `fluff_bus_voltages!`: Process bus voltages to calculate complex voltage forms
- `solution_dictify_buses!`: Transform bus solution data into structured dictionary format
- `solution_dictify_branches!`: Transform branch solution data into structured dictionary format
- `solution_dictify_loads!`: Transform load solution data into structured dictionary format
- `dictify_solution!`: Transform complete power flow solution into structured format

## Base Calculations
- `calc_bases_from_dict`: Calculate electrical base quantities from data dictionary
- `add_vmn_p_q`: Add voltage magnitude, active and reactive power measurements

## Impedance Analysis
- `kron_reduce_impedance`: Perform Kron reduction on 4x4 impedance matrix
- `get_sequence_components`: Calculate sequence components from phase matrix

## Utility Functions
- `show_example`: Display first element of a dictionary
- `show_transformer_math_components`: Extract and display transformer mathematical components

## Result Exploration
- `pf_results`: Get power flow results with summary
- `bus_res`: Get bus results
- `branch_viz`: Visualize branch results

## Engineering Model Exploration
- `eng_report`: Generate engineering model report
- `buses_table`: Display buses table
- `lines_table`: Display lines table
- `loads_table`: Display loads table
- `linecodes_table`: Display linecodes table
- `linecode_table`: Display single linecode details
- `transformers_table`: Display transformers table

## Mathematical Model Exploration
- `math_report`: Generate mathematical model report
- `math_buses_table`: Display mathematical buses table
- `math_branches_table`: Display mathematical branches table
- `math_branch_details`: Display branch details
- `math_loads_table`: Display mathematical loads table
- `math_gen_table`: Display generators table
- `math_gen_details`: Display generator details
- `math_load_details`: Display load details
- `math_bus_details`: Display bus details
- `math_transformers_table`: Display mathematical transformers table
- `math_transformer_details`: Display transformer details

## Network Graph Functions
- `get_graph_node`: Get node properties from network graph
- `get_graph_edge`: Get edge properties from network graph
- `create_network_graph`: Create MetaDiGraph from network data
"""
module PMDUtils

end # module PMDUtils
