using Pkg
Pkg.add(PackageSpec(path=pwd()))

using Pliers
using Pliers.PMDGraph
using Pliers.PMDUtils
using PowerModelsDistribution


# parse dss file

dss_file = joinpath(@__DIR__, "..", "test","data", "trans_example.dss")

eng = parse_file(dss_file)

transform_loops!(eng)

math = transform_data_model(eng, kron_reduce=false,phase_project=false)

plot_network_tree(math, show_node_labels=true, show_edge_labels=true)

math_reduced = deepcopy(math)
reduce_network_buses!(math_reduced)

reduce_network_buses!(math_reduced)


plot_network_tree(math_reduced, show_node_labels=true, show_edge_labels=true)