
# ENG handling and plotting


@testset "ENG handling" begin
    try
        eng = parse_file(pmd_test_models[4][1])
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
        @test true
    catch e
        prinln(e)
        @test false
    end
end

@testset "ENG handling" begin
    for i in eachindex(pmd_test_models)
        eng = parse_file(pmd_test_models[i][1])
        try 
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
            @test true
        catch e
            println(" the error is: ", e)    
            println(" for the network ", pmd_test_models[i][2])
            @test false
        end
    end
end


@testset "Network Graph from ENGINEERING model" begin
    eng = parse_file(pmd_test_models[4][1])
    G, _, d = create_network_graph(eng, !);
    typeof(G)
    @test typeof(G) == MetaDiGraph{Int64, Float64}
    @test typeof(d) == Dict{Symbol, Any}    

    get_graph_node(G, "l1")[:loads][1][:pd_nom]
    @test get_graph_node(G, "l1")[:loads][1][:pd_nom] == [14.0]

    
end