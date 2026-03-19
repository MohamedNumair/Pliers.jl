using LinearAlgebra

@testset "Transformer decomposition reduction" begin
    eng = load_eng_fixture("trans_example.dss")
    PowerModelsDistribution.transform_loops!(eng)
    math = PowerModelsDistribution.transform_data_model(
        eng;
        kron_reduce=false,
        phase_project=false,
    )

    bus_vbase(bus_id) = Float64(get(math["bus"][string(bus_id)], "vbase", get(math["bus"][string(bus_id)], "base_kv", 1.0)))

    grouped_transformers = Dict{String, Vector{String}}()
    for (tx_id, tx) in math["transformer"]
        tx_name = string(get(tx, "name", ""))
        m = match(r"^_virtual_transformer\.(.+)\.(\d+)$", tx_name)
        m === nothing && continue
        group_name = m.captures[1]
        push!(get!(grouped_transformers, group_name, String[]), tx_id)
    end

    selected_group = nothing
    tx_ids = String[]
    vbr_ids = String[]

    for group_name in sort(collect(keys(grouped_transformers)))
        ids = grouped_transformers[group_name]
        branches = [
            br_id
            for (br_id, br) in math["branch"]
            if startswith(string(get(br, "name", "")), "_virtual_branch.transformer.$group_name" * "_")
        ]

        if length(ids) == 2 && !isempty(branches)
            selected_group = group_name
            tx_ids = ids
            vbr_ids = branches
            break
        end
    end

    @test selected_group !== nothing

    if selected_group !== nothing
        virtual_bus_prefix = "_virtual_bus.transformer.$selected_group" * "_"

        ext_vbases = Float64[]
        for tx_id in tx_ids
            tx = math["transformer"][tx_id]
            f_bus = tx["f_bus"]
            t_bus = tx["t_bus"]

            f_name = string(math["bus"][string(f_bus)]["name"])
            t_name = string(math["bus"][string(t_bus)]["name"])

            if startswith(f_name, virtual_bus_prefix)
                push!(ext_vbases, bus_vbase(t_bus))
            elseif startswith(t_name, virtual_bus_prefix)
                push!(ext_vbases, bus_vbase(f_bus))
            end
        end

        @test !isempty(ext_vbases)

        if !isempty(ext_vbases)
            target_vbase = maximum(ext_vbases)

            reduced = deepcopy(math)
            Pliers.PMDUtils.reduce_network_buses_and_transformers!(reduced)

            @test !any(
                startswith(string(get(tx, "name", "")), "_virtual_transformer.$selected_group.")
                for (_, tx) in reduced["transformer"]
            )

            @test !any(
                startswith(string(get(br, "name", "")), "_virtual_branch.transformer.$selected_group" * "_")
                for (_, br) in reduced["branch"]
            )

            eq_ids = [
                br_id
                for (br_id, br) in reduced["branch"]
                if startswith(string(get(br, "name", "")), "_virtual_branch.transformer_reduced.$selected_group")
            ]

            @test length(eq_ids) == 1

            if length(eq_ids) == 1
                eq_branch = reduced["branch"][first(eq_ids)]
                @test isapprox(eq_branch["vbase"], target_vbase; rtol=1e-9, atol=1e-9)
                @test all(isfinite, diag(eq_branch["br_r"]))
                @test all(isfinite, diag(eq_branch["br_x"]))
                @test sum(diag(eq_branch["br_r"])) > 0
                @test sum(diag(eq_branch["br_x"])) > 0
            end
        end
    end
end
