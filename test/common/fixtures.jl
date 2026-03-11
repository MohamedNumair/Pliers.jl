const DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))
const DSS_FIXTURES = ("ieee-33-bus.dss", "trans_example.dss")

fixture_path(name::AbstractString) = normpath(joinpath(DATA_DIR, name))

function load_eng_fixture(name::AbstractString)
    path = fixture_path(name)
    isfile(path) || error("Missing test fixture: $path")
    return PowerModelsDistribution.parse_file(path)
end

function capture_stdout(f::Function)
    return Suppressor.@capture_out begin
        f()
    end
end
