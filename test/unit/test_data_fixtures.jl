@testset "Local fixture files" begin
    for fixture in DSS_FIXTURES
        @test isfile(fixture_path(fixture))
    end
end
