@testset "PMDUtils placeholders" begin
    @testset "Reporting APIs" begin
        eng = load_eng_fixture("trans_example.dss")

        eng_report_out = capture_stdout() do
            Pliers.PMDUtils.eng_report(eng)
        end
        @test occursin("This network has", eng_report_out)

        eng_report_detailed_out = capture_stdout() do
            Pliers.PMDUtils.eng_report(eng; detailed=true)
        end
        @test occursin("Buses Table", eng_report_detailed_out)
        @test occursin("Lines Table", eng_report_detailed_out)
        @test occursin("Loads Table", eng_report_detailed_out)

        math_eng = load_eng_fixture("ieee-33-bus.dss")
        math = PowerModelsDistribution.transform_data_model(math_eng; kron_reduce=false, phase_project=false)
        math_report_out = capture_stdout() do
            Pliers.PMDUtils.math_report(math)
        end
        @test occursin("This network has", math_report_out)

        math_report_detailed_out = capture_stdout() do
            Pliers.PMDUtils.math_report(math; detailed=true)
        end
        @test occursin("Buses Table", math_report_detailed_out)
        @test occursin("Branches Table", math_report_detailed_out)
        @test occursin("Loads Table", math_report_detailed_out)
    end

    @testset "Filtering table APIs" begin
        eng = load_eng_fixture("trans_example.dss")

        sample_bus_id = first(keys(eng["bus"]))
        matched_buses = Pliers.PMDUtils.buses_table(deepcopy(eng), bus -> bus["bus_id"] == sample_bus_id)
        @test sample_bus_id in matched_buses

        connected_lines = Pliers.PMDUtils.lines_table(deepcopy(eng), line -> line["f_bus"] == sample_bus_id || line["t_bus"] == sample_bus_id)
        @test !isempty(connected_lines)

        matched_loads = Pliers.PMDUtils.loads_table(deepcopy(eng), load -> get(load, "connections", Int[]) == [1, 4])
        @test !isempty(matched_loads)

        linecodes_out = capture_stdout() do
            Pliers.PMDUtils.linecodes_table(eng)
        end
        @test occursin("Linecodes Table", linecodes_out)

        math_eng = load_eng_fixture("ieee-33-bus.dss")
        math = PowerModelsDistribution.transform_data_model(math_eng; kron_reduce=false, phase_project=false)
        sample_f_bus = first(values(math["branch"]))["f_bus"]
        matched_branches = Pliers.PMDUtils.math_branches_table(math, branch -> branch["f_bus"] == sample_f_bus)
        @test !isempty(matched_branches)

        matched_math_loads = Pliers.PMDUtils.math_loads_table(math, load -> !isempty(load["connections"]))
        @test !isempty(matched_math_loads)

        # TODO: `math_branch_details` currently errors with PrettyTables header keyword handling.
        @test_skip false
        # TODO: `math_load_details` currently errors with PrettyTables header keyword handling.
        @test_skip false
    end

    @testset "Network reduction APIs" begin
        # TODO: cover reduce_network_intermediate_buses!, reduce_empty_leaf_buses!, reduce_network_buses!
        @test_skip false
    end

    @testset "Direction and topology fixes" begin
        # TODO: cover fix_eng_directions! and related edge-orientation behavior.
        @test_skip false
    end

    @testset "Solution shaping helpers" begin
        # TODO: cover fluff_bus_voltages! and dictify_solution! family on realistic PF outputs.
        @test_skip false
    end

    @testset "Impedance helpers" begin
        # TODO: cover kron_reduce_impedance and get_sequence_components with deterministic fixtures.
        @test_skip false
    end
end
