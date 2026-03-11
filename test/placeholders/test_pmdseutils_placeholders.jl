@testset "PMDSEUtils placeholders" begin
    @testset "Residual summary tables" begin
        # TODO: validate df_meas_res and viz_residuals output schema for a deterministic SE case.
        @test_skip false
    end

    @testset "Measurement augmentation" begin
        # TODO: test add_pd_qd_vmn! and write_sm_measurements side effects on math data.
        @test_skip false
    end
end
