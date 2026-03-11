using Graphs

@testset "Network graph creation" begin
    eng = load_eng_fixture("ieee-33-bus.dss")
    eng_graph, eng_layout, eng_sym = Pliers.PMDGraph.create_network_graph(eng, Pliers.PMDGraph.smart_layout)

    @test Graphs.nv(eng_graph) > 0
    @test Graphs.ne(eng_graph) > 0
    @test isa(eng_sym, Dict{Symbol,Any})
    @test eng_layout isa Function || eng_layout === Pliers.PMDGraph.smart_layout

    math = PowerModelsDistribution.transform_data_model(eng; kron_reduce=false, phase_project=false)
    math_graph, math_layout, math_sym = Pliers.PMDGraph.create_network_graph(math, Pliers.PMDGraph.smart_layout)

    @test Graphs.nv(math_graph) > 0
    @test Graphs.ne(math_graph) > 0
    @test isa(math_sym, Dict{Symbol,Any})
    @test math_layout isa Function || math_layout === Pliers.PMDGraph.smart_layout
end
