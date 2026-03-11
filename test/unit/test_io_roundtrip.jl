@testset "Serialization IO" begin
    payload = Dict("name" => "fixture", "values" => [1, 2, 3])

    mktemp() do path, io
        close(io)
        Pliers.save_data(payload, path)
        loaded = Pliers.read_data(path)
        @test loaded == payload
    end
end
