"""
Example: plot_network_by_voltage with trans_example.dss
=========================================================
This script loads the transformer test network (two voltage zones: 11 kV and
0.4 kV), calls `plot_network_by_voltage`, and prints a detailed diagnostic
report so you can see exactly which vbase value was resolved for every bus.

Run with:
    julia --project=. Examples/plot_voltage_example.jl

Enable debug logging to see the internal resolution steps:
    JULIA_DEBUG=all julia --project=. Examples/plot_voltage_example.jl
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Logging
using PowerModelsDistribution
using Pliers
import Pliers.PMDGraph:
    _extract_vbase_map_eng,
    _extract_vbase_map_math,
    _format_voltage_label

# ---------------------------------------------------------------------------
# 1.  Enable debug logging so every internal step is visible
# ---------------------------------------------------------------------------
global_logger(ConsoleLogger(stderr, Logging.Debug))

# ---------------------------------------------------------------------------
# 2.  Load the network
# ---------------------------------------------------------------------------
dss_file = joinpath(@__DIR__, "..", "test", "data", "trans_example.dss")
@info "Loading DSS file: $dss_file"
eng = PowerModelsDistribution.parse_file(dss_file)

@info "data_model = $(eng["data_model"])"
@info "Buses in engineering model: $(collect(keys(eng["bus"])))"

# ---------------------------------------------------------------------------
# 3.  Inspect raw engineering bus data
# ---------------------------------------------------------------------------
println("\n=== Raw engineering bus keys ===")
for (bus_id, bus_data) in eng["bus"]
    println("  bus '$bus_id'  → keys: $(sort(collect(keys(bus_data))))")
    for (k, v) in bus_data
        println("      $k : $v")
    end
end

# ---------------------------------------------------------------------------
# 4.  Inspect raw transformer data (this is where vm_nom lives)
# ---------------------------------------------------------------------------
println("\n=== Raw transformer data ===")
if haskey(eng, "transformer")
    for (tid, trans) in eng["transformer"]
        println("  transformer '$tid'")
        for (k, v) in trans
            println("      $k : $v")
        end
    end
else
    println("  (no transformers found)")
end

# ---------------------------------------------------------------------------
# 5.  Inspect raw voltage_source data
# ---------------------------------------------------------------------------
println("\n=== Raw voltage_source data ===")
for (vsid, vs) in eng["voltage_source"]
    println("  voltage_source '$vsid'")
    for (k, v) in vs
        println("      $k : $v")
    end
end

# ---------------------------------------------------------------------------
# 6.  Run vbase extraction manually and print results
# ---------------------------------------------------------------------------
using Pliers: convert_keys_to_symbols
eng_sym = convert_keys_to_symbols(eng)

println("\n=== _extract_vbase_map_eng output ===")
vbase_map = _extract_vbase_map_eng(eng_sym)
if isempty(vbase_map)
    println("  WARNING: vbase_map is EMPTY — all buses will show as 'Unknown'")
else
    for (bus_id, kv) in sort(collect(vbase_map), by=first∘string)
        println("  bus :$bus_id → $(_format_voltage_label(kv))")
    end
end

# ---------------------------------------------------------------------------
# 7.  Optionally: run on the math model too
# ---------------------------------------------------------------------------
println("\n=== Math model vbase (transform_data_model) ===")
math = PowerModelsDistribution.transform_data_model(eng)
math_sym = convert_keys_to_symbols(math)
vbase_map_math = _extract_vbase_map_math(math_sym)
if isempty(vbase_map_math)
    println("  WARNING: math vbase_map is EMPTY")
else
    for (bus_id, kv) in sort(collect(vbase_map_math), by=first∘string)
        println("  math bus :$bus_id → $(_format_voltage_label(kv))")
    end
end

# ---------------------------------------------------------------------------
# 8.  Call the actual plotting function (engineering model)
# ---------------------------------------------------------------------------
println("\n=== Calling plot_network_by_voltage (engineering model) ===")
using CairoMakie   # use CairoMakie for non-interactive / file output
CairoMakie.activate!()

fig_eng = plot_network_by_voltage(
    eng;
    makie_backend=CairoMakie,
    figure_size=(1200, 800),
    show_node_labels=true,
)
output_eng = joinpath(@__DIR__, "voltage_plot_eng.png")
save(output_eng, fig_eng)
@info "Saved engineering model plot → $output_eng"

# ---------------------------------------------------------------------------
# 9.  Call the actual plotting function (math model)
# ---------------------------------------------------------------------------
println("\n=== Calling plot_network_by_voltage (math model) ===")
fig_math = plot_network_by_voltage(
    math;
    makie_backend=CairoMakie,
    figure_size=(1200, 800),
    show_node_labels=true,
)
output_math = joinpath(@__DIR__, "voltage_plot_math.png")
save(output_math, fig_math)
@info "Saved math model plot → $output_math"

println("\nDone. Check $(dirname(output_eng)) for the PNG files.")
