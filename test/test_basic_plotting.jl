using Test

@testset "Basic plotting smoke" begin
    eng = load_eng_fixture("trans_example.dss")

    f_tree, ax_tree, p_tree = Pliers.PMDGraph.plot_network_tree(
        eng;
        makie_backend=Pliers.CairoMakie,
        show_node_labels=true,
        show_edge_labels=true,
    )
    @test !isnothing(f_tree)
    @test !isnothing(ax_tree)
    @test !isnothing(p_tree)

    out_png = tempname() * ".png"
    Pliers.CairoMakie.save(out_png, f_tree)
    @test isfile(out_png)
    @test filesize(out_png) > 0

    f_coords, ax_coords, p_coords = Pliers.PMDGraph.plot_network_coords(
        eng;
        makie_backend=Pliers.CairoMakie,
        show_node_labels=false,
        show_edge_labels=false,
        show_load_labels=true,
    )
    @test !isnothing(f_coords)
    @test !isnothing(ax_coords)
    @test !isnothing(p_coords)
end
