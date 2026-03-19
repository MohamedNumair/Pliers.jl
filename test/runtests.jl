using Test

using Pliers
using Pliers.PMDUtils
using Pliers.PMDGraph
using PowerModelsDistribution
using Suppressor

include("common/fixtures.jl")

@testset "Pliers.jl" begin
    @testset "Fixture Integrity" begin
        include("unit/test_data_fixtures.jl")
    end

    @testset "Core Smoke Tests" begin
        include("unit/test_io_roundtrip.jl")
        include("unit/test_parse_and_transform.jl")
        include("unit/test_graph_smoke.jl")
        include("unit/test_math_reduction_transformers.jl")
    end

    @testset "Coverage Placeholders" begin
        include("test_basic_plotting.jl")
        include("placeholders/test_pmdutils_placeholders.jl")
        include("placeholders/test_pmdgraph_placeholders.jl")
        include("placeholders/test_pmdseutils_placeholders.jl")
        include("placeholders/test_io_placeholders.jl")
    end
end




