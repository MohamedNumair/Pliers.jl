@testset "Dataset loader placeholders" begin
    @testset "ENWL dataset loaders" begin
        data_1ph_path = joinpath(Pliers.BASE_DIR, "test", "data", "ENWL1ph", "Ntw_1_Fdr_1.jld2")
        @test isfile(data_1ph_path)

        eng_1ph = Pliers.load_enwl_model(1, 1, "1ph")
        @test isa(eng_1ph, Dict{String,Any})
        @test haskey(eng_1ph, "bus")
        @test !isempty(eng_1ph["bus"])
    end

    @testset "Spanish dataset loaders" begin
        spanish_dataset_path = joinpath(Pliers.BASE_DIR, "test", "data", "Spanish_Network.jld2")
        @test isfile(spanish_dataset_path)

        spanish_dataset = Pliers.load_spanish_dataset()
        @test isa(spanish_dataset, Dict)

        network_strings = Pliers.spanish_network_strings()
        @test !isempty(network_strings)

        network_data = Pliers.load_spanish_network(first(network_strings))
        @test isa(network_data, Dict) || isa(network_data, Vector)
    end

    @testset "Path extraction helpers" begin
        ntw, fdr = Pliers.extract_network_and_feeder(raw"C:\Data\Network_17\Feeder_4\Master.dss")
        @test ntw == 17
        @test fdr == 4

        ntw_jld, fdr_jld = Pliers.extract_ntw_fdr_jld(raw"C:\Data\MENWL\Ntw_18_Fdr_9.jld2")
        @test ntw_jld == 18
        @test fdr_jld == 9
    end
end
