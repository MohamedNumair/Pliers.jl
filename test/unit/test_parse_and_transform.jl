@testset "DSS parsing" begin
    for fixture in DSS_FIXTURES
        eng = load_eng_fixture(fixture)
        @test isa(eng, Dict{String,Any})
        @test haskey(eng, "data_model")
        @test eng["data_model"] == PowerModelsDistribution.ENGINEERING
        @test haskey(eng, "bus")
        @test !isempty(eng["bus"])
    end
end

@testset "ENG to MATH transform" begin
    eng = load_eng_fixture("ieee-33-bus.dss")
    math = PowerModelsDistribution.transform_data_model(eng; kron_reduce=false, phase_project=false)

    @test isa(math, Dict{String,Any})
    @test haskey(math, "data_model")
    @test math["data_model"] == PowerModelsDistribution.MATHEMATICAL
    @test haskey(math, "bus")
    @test haskey(math, "branch")
    @test !isempty(math["bus"])
    @test !isempty(math["branch"])
end
