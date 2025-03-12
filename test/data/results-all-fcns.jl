using Pliers
using Pliers.PowerModelsDistribution

using Pliers.Proj
using Pliers.CSV
using Pliers.DataFrames
using Pliers.FileIO


# Loading Test files
# Load the model
pmd_path = joinpath(dirname(pathof(PowerModelsDistribution)), "..")

file_options= [
    joinpath(pmd_path, "test/data/opendss/case3_balanced.dss") => "case3_balanced",
    joinpath(pmd_path, "test/data/opendss/case3_unbalanced.dss") => "case3_unbalanced",
    joinpath(pmd_path, "test/data/opendss/case3_balanced_battery.dss") => "case3_balanced_battery",
    joinpath(pmd_path, "test/data/opendss/case5_phase_drop.dss") => "case5_phase_drop",
    joinpath(pmd_path, "test/data/opendss/ut_trans_2w_yy_oltc.dss") => "ut_trans_2w_yy_oltc",
    joinpath(pmd_path, "test/data/opendss/case3_balanced_battery.dss") => "case3_balanced_battery",
]
en_cases_files = readdir(joinpath(pmd_path, "test/data/en_validation_case_data"); join=true)

joinpath(joinpath(dirname(pathof(PowerModelsDistribution)), ".."), "test/data/opendss/case3_unbalanced.dss")
#enwl_en_file = readdir("C:/Users/mnumair/OneDrive - KU Leuven/Research/Distribution Networks/Four-wire_low_voltage_power_network_dataset-QEzDvqEq-/data/Four-wire"; join=true)

# for file in en_cases_files
#     println("$(file)")
#     eng= parse_file(file)

#     eng_report(eng)
#     plot_network_tree(eng, makie_backend=WGLMakie)
# end

eng= parse_file(file_options[4][1])
eng= parse_file(en_cases_files[1])
# or 

enwl_1_1 = load_enwl_model(1,1)
enwl_1_2 = load_enwl_model(1,2)
enwl_1_3 = load_enwl_model(1,3)
enwl_1_4 = load_enwl_model(1,4)
eng = enwl_1_1

math = transform_data_model(eng, kron_reduce=false)


math_kron = transform_data_model(eng, kron_reduce=true)

# ENG handling and plotting

eng_report(eng)
eng_report(eng, detailed= true)
buses_table(eng)
buses_table(eng, bus -> bus["bus_id"] =="sourcebus")
lines_table(eng)
lines_table(eng, line -> line["f_bus"] == "sourcebus")
lines_table(eng, line -> line["length"] > 0.75)
loads_table(eng)
loads_table(eng, load -> load["pd_nom"] > [0.33])
loads_table(eng, load -> load["connections"] == [1, 4])
loads_table(eng, load -> load["connections"] == [2, 4])
loads_table(eng, load -> load["connections"] == [3, 4])
linecodes_table(eng)


G, _, eng_sym = create_network_graph(eng, !);
G, _, _ = create_network_graph(enwl_1_1, !);


for (_,edge) in G.eprops
    display(edge)
end


get_graph_node(G, "l1")[:loads]
get_graph_node(G, "sourcebus", "voltage_sources")

G[1, :bus_id]


get_graph_edge(G, "line2")
get_graph_edge(G, "line1")[:linecodes][:rs]
get_graph_edge(G, "line1","linecodes")[:rs]
get_graph_edge(G, "line1","f_bus")


ne(G)


findfirst(v -> !isempty(v[:loads]), G.vprops)

G.vprops[4][:loads][2][:connections]


for (n, node) in G.vprops
    display(node)

end

for (e, edge) in G.eprops
    display(edge[:line_id])
end


plot_network_tree(enwl_1_1,
                     show_node_labels = false,
                     show_edge_labels = false,
                     makie_backend=Pliers.WGLMakie
                )

plot_network_coords(eng, show_edge_labels = true)
plot_network_coords!(eng2, show_edge_labels = false)


f, ax, p = plot_network_coords(enwl_fix_1_1,
                               show_edge_labels = true,
                               #edge_color = :black,
                               #node_color = :red,
                               #node_size = 1,
                               makie_backend=GLMakie,
                               size = (80, 5000)
                               
                            )
                            hidespines!(ax)
                            hidedecorations!(ax)
                            f
                            # Legend(f[1,1], [LineElement(color = :red, linestyle = nothing),
                            #                 LineElement(color = :green, linestyle = nothing),
                            #                 LineElement(color = :blue, linestyle = nothing),
                            #                 LineElement(color = :purple, linestyle = nothing),],
                            #                 ["phase cable (a-n)", "phase cable (b-n)", "phase cable (c-n)", "Three-phase cable (abc-n)" 
                            
                            #                 ],)
                            #                 #patchsize = (35, 35), rowgap = 10)
                            

                            save("fixed_files_lines.svg", f)
                            





plot_network_map(
                    enwl_1_1,                     
                    makie_backend = Pliers.GLMakie,
                   # node_size = 1,
                    #edge_color = :red,
                    #node_color = :red,
                    show_node_labels = false,
                    show_edge_labels = false,
                    edge_labels_type = :line_id,
                    figure_size = (5000,5000),
                    tiles_provider =  Pliers.TileProviders.Google(:satelite), # :roadmap, :satelite, :terrain, :hybrid

                    
                )
plot_network_coords!(enwl_1_2,
                     show_edge_labels = false,
                     makie_backend=Pliers.GLMakie,
                     edge_color = :blue,
                     #node_color = :blue,
                    # node_size = 1,
                    )
plot_network_coords!(enwl_1_3,
                     show_edge_labels = false,
                     makie_backend=Pliers.GLMakie,
                     edge_color = :purple,
                     #node_color = :purple,
                    # node_size = 1,
                    )
plot_network_coords!(enwl_1_4,
                     show_edge_labels = false,
                     makie_backend=Pliers.GLMakie,
                     edge_color = :green,
                     
                     # node_color = :green,
                     #node_size = 1,
                    )

hidedecorations!(ax); 
f
save("ENWL_ntw1.pdf", f)

# MATH handling

math = transform_data_model(eng)

math = transform_data_model(eng, kron_reduce=false, phase_project=false)
add_start_vrvi!(math)

math_report(math)   
math_report(math, detailed = true)   
math_buses_table(math)
math_buses_table(math, bus -> bus["bus_type"] == 3)
math_bus_details(math, ["5"])
math_branches_table(math)
idx = math_branches_table(math, branch -> branch["f_bus"] == 2); 
math_branch_details(math, idx)
math_branch_details(math, ["5"])
idx = math_loads_table(math, load -> load["connections"] == [3])
math_load_details(math, idx)  


math_gen_table(math, gen -> gen["cost"][1] > 50) 

math_gen_details(math, ["1"])


plot_network_tree(math,  show_node_labels = true, show_edge_labels=true, makie_backend=GLMakie)                             
plot_network_tree(eng, show_node_labels = true,show_edge_labels = true, makie_backend=GLMakie)                             


G_math, _, _ = create_network_graph(math, !);
    
get_graph_edge(G_math, "1")


# Compare both CSIRO 4-wire and my 4-wire models


csiro_11 = parse_file("C:/Users/mnumair/OneDrive - KU Leuven/PhD Files/2- Repositories and Models/Four-wire_low_voltage_power_network_dataset-QEzDvqEq-/data/Four-wire/network_1/Feeder_1/Master.dss")
transform_loops!(csiro_11)
loads_table(csiro_11)
math_cisro_11 = transform_data_model(csiro_11, kron_reduce=false, phase_project=false)  
add_start_vrvi!(math_cisro_11)
PF_RES_math_cisro_11 = PowerModelsDistribution.compute_mc_pf(math_cisro_11; explicit_neutral=true, max_iter=15)
pf_sol_cs = PF_RES_math_cisro_11["solution"]


fixed4_11 = parse_file("C:/Users/mnumair/OneDrive - KU Leuven/PhD Files/2- Repositories and Models/ENWL original dataset/lv-network-models-2/enwl_og_networks_corrected_linePhase_fourWire/network_1/Feeder_1/Master.dss") 
loads_table(fixed4_11)
plot_network_coords(fixed4_11, show_node_labels = true, makie_backend=GLMakie)
transform_loops!(fixed4_11)
math_fixed4_11 = transform_data_model(fixed4_11, kron_reduce=false, phase_project=false)

add_start_vrvi!(math_fixed4_11)
PF_RES_math_fixed4_11 = PowerModelsDistribution.compute_mc_pf(math_fixed4_11; explicit_neutral=true, max_iter=15)
pf_sol_fx = PF_RES_math_fixed4_11["solution"]


pf_sol_cs["bus"]["59"]
pf_sol_fx["bus"]["59"]

math_branches_table(math_fixed4_11, branch -> branch["source_id"] == "line.line58")
math_branches_table(math_fixed4_11, branch -> branch["source_id"] == "line.line64")
math_branches_table(math_fixed4_11, branch -> branch["source_id"] == "line.line65")

pf_sol_cs["branch"]["429"] # line 58
pf_sol_cs["branch"]["62"] # line 64
pf_sol_cs["branch"]["542"] # line 65

pf_sol_fx["branch"]["429"] # line 58
pf_sol_fx["branch"]["62"] # line 64
pf_sol_fx["branch"]["542"] # line 65


math_cisro_11["branch"]["542"]["br_r"] # line 58
math_cisro_11["branch"]["542"]["br_x"] # line 58

math_fixed4_11["branch"]["542"]["br_r"]
math_fixed4_11["branch"]["542"]["br_x"]


csiro_11["linecode"][csiro_11["line"]["line65"]["linecode"]]
fixed4_11["linecode"][fixed4_11["line"]["line65"]["linecode"]]


csiro_11["linecode"][csiro_11["line"]["line65"]["linecode"]]["rs"]
fixed4_11["linecode"][fixed4_11["line"]["line65"]["linecode"]]["rs"]

csiro_11["linecode"][csiro_11["line"]["line65"]["linecode"]]["xs"]
fixed4_11["linecode"][fixed4_11["line"]["line65"]["linecode"]]["xs"]






fixed4_11["linecode"]["lc1"]["rs"]

fixed4_11

math_cisro_11["branch"]["542"]["br_r"]
math_cisro_11["branch"]["542"]["br_x"]

math_fixed4_11["branch"]["542"]["br_r"]
math_fixed4_11["branch"]["542"]["br_x"]



Ic_csiro = pf_sol_cs["branch"]["542"]["cr"] .+ im*pf_sol_cs["branch"]["542"]["cr"] # line 65

abs.(Ic_csiro)  


pf_sol_fx["branch"]["542"] # line 65

ic_fx = pf_sol_fx["branch"]["542"]["cr"] .+ im*pf_sol_fx["branch"]["542"]["cr"] # line 65
abs.(ic_fx)

##############################################################

#Line Geometry

lG = parse_file("C:/Users/mnumair/OneDrive - KU Leuven/PhD Agenda/2024/24-12 December/Line Geometry/LineExampleGeometry.dss")

buses_table(lG)
lines_table(lG)

lG["line"]["lineexample"]["rs"]

testGeometry = parse_file("C:/Users/mnumair/OneDrive - KU Leuven/PhD Agenda/2024/24-12 December/Line Geometry/LineExampleGeometry.dss")

buses_table(testGeometry)
lines_table(testGeometry)

testGeometry["line"]["lineexample"]["rs"]
testGeometry["line"]["lineexample"]["xs"]

math = transform_data_model(testGeometry, kron_reduce=false, phase_project=false)

ieeeG = parse_file("C:/Users/mnumair/OneDrive - KU Leuven/PhD Agenda/2024/24-12 December/Line Geometry/IEEEGeometries.dss")



ieeeG["bus"]


spain = parse_file("C:/Users/mnumair/OneDrive - KU Leuven/PhD Agenda/2024/24-11 November/TSK-478 mappingEULV/Mapping Spanish/Master.dss")
transform_loops!(spain)
remove_all_bounds!(spain)
spain["settings"]["base_frequency"] = 50
Pliers.add_feeders_cktbks!(spain)



feeder = Pliers.isolate_feeder(spain, 160)




for fno in 1:160
    Pliers.write_dss_file(spain, fno)
end







for fno in 1:160
    feeder = parse_file("C:/Users/mnumair/OneDrive - KU Leuven/PhD Agenda/2024/24-11 November/TSK-478 mappingEULV/Mapping Spanish/spanish_feeders/Feeder_$fno.dss")
    transform_loops!(feeder)
    remove_all_bounds!(feeder)
    Pliers.FileIO.save("C:/Users/mnumair/OneDrive - KU Leuven/PhD Agenda/2024/24-11 November/TSK-478 mappingEULV/Mapping Spanish/spanish_feeders_eng/Feeder_$fno.jld2", feeder)
end

directory = "C:/Users/mnumair/OneDrive - KU Leuven/PhD Agenda/2024/24-11 November/TSK-478 mappingEULV/Mapping Spanish/spanish_feeders"


file_names = []
files_path = []

for (root, dirs, file) in walkdir(directory)
    for f in file
        if occursin(".dss", f)
            push!(file_names, f)
            push!(files_path, joinpath(root, f))
        end
    end
end

files_path
file_names

for (file_name, file_path) in zip(file_names, files_path)
    println(file_name)
    println(file_path)
    feeder= parse_file(file_path)
    transform_loops!(feeder)
    remove_all_bounds!(feeder)
    file_name_first= split(file_name, ".")[1]
    Pliers.FileIO.save("C:/Users/mnumair/OneDrive - KU Leuven/PhD Agenda/2024/24-11 November/TSK-478 mappingEULV/Mapping Spanish/spanish_feeders_eng/$file_name_first.jld2", feeder)

end


transform_loops!(Feeder_1_spanish)


plot_network_map(Feeder_1_spanish)

math = transform_data_model(Feeder_1_spanish, kron_reduce=false, phase_project=false)


feeder = Pliers.load_spanish_feeder(1)

feeder[1]

feeders = Pliers.load_spanish_network("td400291")


Pliers.network_strings()[1]


Network_number = 29
feeders = Pliers.load_spanish_network(Pliers.network_strings()[Network_number])


plot_network_map(feeders[1])

for feeder in feeders 
    plot_network_coords!(feeder)
end



f1 = Pliers.load_spanish_feeder(1)

plot_network_map(f1)

for i in 2:160
    f = Pliers.load_spanish_feeder(i)
    plot_network_coords!(f)
end
