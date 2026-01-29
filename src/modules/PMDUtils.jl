"""
    PMDUtils

Internal sub-module providing utility functions for PowerModelsDistribution (PMD) workflows.

This module re-exports functions from the main Pliers module for:
- Processing power flow solutions (voltage fluffing, dictification)
- Impedance calculations (Kron reduction, sequence components)
- Network data manipulation
- Result analysis and transformation
- Engineering and mathematical model exploration

See the main Pliers module for function documentation.
"""
module PMDUtils

using ..Pliers

# Re-export PMD utility functions from parent module
export fluff_bus_voltages!
export solution_dictify_buses!, solution_dictify_branches!, solution_dictify_loads!, dictify_solution!
export calc_bases_from_dict, add_vmn_p_q
export kron_reduce_impedance, get_sequence_components
export show_example, show_transformer_math_components

# Re-export result explorer functions
export pf_results, bus_res, branch_viz

# Re-export eng explorer functions
export eng_report, buses_table, lines_table, loads_table, linecodes_table, linecode_table, transformers_table

# Re-export math explorer functions
export math_report, math_buses_table, math_branches_table, math_branch_details
export math_loads_table, math_gen_table, math_gen_details, math_load_details
export math_bus_details, math_transformers_table, math_transformer_details

# Re-export network graph functions
export get_graph_node, get_graph_edge, create_network_graph

end # module PMDUtils
