@testset "PMDGraph placeholders" begin
    @testset "Graph semantics" begin
        eng = load_eng_fixture("ieee-33-bus.dss")
        eng_graph, _, _ = Pliers.PMDGraph.create_network_graph(eng, Pliers.PMDGraph.smart_layout)

        bus_id = first(keys(eng["bus"]))
        node_props = Pliers.PMDGraph.get_graph_node(eng_graph, bus_id)
        @test get(node_props, :bus_id, nothing) == Symbol(bus_id)
        @test Pliers.PMDGraph.get_graph_node(eng_graph, bus_id, "bus_id") == Symbol(bus_id)

        line_id = first(keys(eng["line"]))
        edge_props = Pliers.PMDGraph.get_graph_edge(eng_graph, line_id)
        @test get(edge_props, :line_id, nothing) == Symbol(line_id)
        @test Pliers.PMDGraph.get_graph_edge(eng_graph, line_id, "line_id") == Symbol(line_id)

        math = PowerModelsDistribution.transform_data_model(eng; kron_reduce=false, phase_project=false)
        math_graph, _, _ = Pliers.PMDGraph.create_network_graph(math, Pliers.PMDGraph.smart_layout)

        branch_id = first(keys(math["branch"]))
        math_edge_props = Pliers.PMDGraph.get_graph_edge(math_graph, branch_id)
        @test get(math_edge_props, :branch_id, nothing) == Symbol(branch_id)
    end

    @testset "Tree plotting smoke tests" begin
        eng = load_eng_fixture("ieee-33-bus.dss")
        f, ax, p = Pliers.PMDGraph.plot_network_tree(
            eng;
            makie_backend=Pliers.CairoMakie,
            show_node_labels=true,
            show_edge_labels=true,
        )
        @test !isnothing(f)
        @test !isnothing(ax)
        @test !isnothing(p)

        out_png = tempname() * ".png"
        Pliers.CairoMakie.save(out_png, f)
        @test isfile(out_png)
        @test filesize(out_png) > 0
    end

    @testset "Coordinate plotting smoke tests" begin
        eng = load_eng_fixture("ieee-33-bus.dss")
        f, ax, p = Pliers.PMDGraph.plot_network_coords(
            eng;
            makie_backend=Pliers.CairoMakie,
            show_node_labels=false,
            show_edge_labels=false,
            show_load_labels=true,
        )
        @test !isnothing(f)
        @test !isnothing(ax)
        @test !isnothing(p)
    end

    @testset "smart_layout behavior" begin
        # TODO: ensure smart_layout remains stable on cyclic and disconnected cases.
        @test_skip false
    end
end
