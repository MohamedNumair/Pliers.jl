# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Activate environment
julia --project=.

# Run all tests
julia --project=. test/runtests.jl
# or in Julia REPL: pkg> test

# Run a single test file
julia --project=. -e 'include("test/unit/test_parse_and_transform.jl")'

# Build documentation
julia --project=docs docs/make.jl
```

## Architecture

Pliers.jl is a toolkit for power distribution network analysis, built on top of `PowerModelsDistribution.jl` (PMD) and `PowerModelsDistributionStateEstimation.jl` (PMDSE). All network data is represented as `Dict{String, Any}`.

**Three main modules** under `src/modules/`, re-exported from `src/Pliers.jl`:

- **`PMDUtils`** (~3,700 lines): Network reporting (buses, lines, loads, transformers as `DataFrame`/`PrettyTable`), result processing (power flow solutions, voltage profiles), impedance calculations (Kron reduction, sequence components), and solution dictification.
- **`PMDGraph`** (~2,300 lines): Graph creation (`create_network_graph()`), tree layout (`smart_layout()`), and visualization — `plot_network_tree()`, `plot_network_coords()`, `plot_network_map()`, and `bus_phasor()`.
- **`PMDSEUtils`** (~600 lines): State estimation result visualization and measurement residual analysis.

**Supporting files:**
- `src/core/export.jl`: Centralized exports for all public functions
- `src/core/styles.jl`: Makie plot themes for journal-quality output
- `src/core/utils.jl`: General helpers (key extraction, vector comparison)
- `src/io/`: Network loading/saving for ENWL (UK) and Spanish feeder datasets

## Multi-Root Workspace

This repository is one of four roots in the active workspace:

| Root | Path | Purpose |
|------|------|---------|
| **Pliers** | `c:\Users\mnumair\.julia\dev\Pliers\` | This package — analysis, reporting, visualization |
| **TSK-844** | `c:\onedrive\PhD Agenda\2026\03-26 March\TSK-844 Timerseries MVLV transformers\` | Active research — MV/LV transformer monitoring via timeseries DSSE |
| **PMD** | `c:\Users\mnumair\.julia\dev\PowerModelsDistribution\` | Upstream — distribution network modelling (v0.16.0) |
| **PMDSE** | `c:\Users\mnumair\.julia\dev\PowerModelsDistributionStateEstimation\` | Upstream — DSSE (v0.8.0, branch: `feat/transformer_handling`) |

**Quick reference:**

| Short name | Julia usage | When |
|------------|-------------|------|
| PMD | `import PowerModelsDistribution as PMD` | Upstream modelling functions |
| PMDSE | `import PowerModelsDistributionStateEstimation` | SE functions |
| ENG model | `data = parse_file(dss_path)` | String-keyed engineering data |
| MATH model | `math = transform_data_model(eng)` | Integer-keyed mathematical data |
| MN data | `mn_data = make_multinetwork(data)` | Timeseries / multi-network structures |

## PMD — Data Model Pipeline & Source Structure

```
OpenDSS (.dss) → rawdss2dss → dss2eng → ENG model (Dict{String,Any}, String keys)
                                              ↓ transform_data_model()
                                   MATH model (Dict{String,Any}, Int keys)
                                              ↓ math2eng
                                   Solution mapped back to ENG
```

**Key source files** (`src/data_model/transformations/`):

- `rawdss2dss.jl` — Raw DSS → structured DSS
- `dss2eng.jl` — DSS → Engineering model
- `eng2math.jl` — Engineering → Mathematical model
- `math2eng.jl` — Mathematical → Engineering (solution mapping)
- `reduce.jl` — Network reduction

**Formulations** (`src/form/`): ACP, ACR, IVR, LinDist3Flow, and others.

### Key PMD Functions

- `parse_file(path)` — parse OpenDSS → ENG model
- `transform_data_model(eng)` — ENG → MATH
- `make_multinetwork(data)` — create timeseries multi-network structure
- `solve_mc_opf(data, form, solver)` — multi-conductor OPF
- `solve_mc_pf(data, form, solver)` — multi-conductor power flow

### ENG vs MATH Differences

| Aspect | ENG | MATH |
|--------|-----|------|
| Indexing | String IDs | Integer indices |
| Lines | `data["line"]` | `data["branch"]` |
| Transformers | `data["transformer"]` | Decomposed: 2 ideal transformers + 3 virtual branches |
| Virtual elements | None | `_virtual_bus`, `_virtual_branch` added |

Always verify which model you are working with — indexing and component names differ.

## PMDSE — Source Structure & Key Functions

**Source layout:**

```
src/
├── bad_data/
│   ├── chi_squares_test.jl                  # Chi-squared hypothesis test
│   └── largest_normalized_residuals.jl      # LNR bad data identification
├── core/
│   ├── constraint.jl                        # SE-specific constraints
│   ├── objective.jl                         # WLS, WLAV, LAV objectives
│   ├── variable.jl                          # SE variables (measurement residuals)
│   ├── measurement_conversion.jl            # Convert measurements between component types
│   └── start_values.jl                      # Initial value computation
├── form/
│   ├── adapted_pmd_constraints.jl           # Modified PMD constraints for SE
│   ├── reduced_ac.jl                        # Reduced AC formulations
│   └── reduced_ivr.jl                       # Reduced IVR formulations
├── io/
│   ├── measurements.jl                      # Measurement parsing and configuration
│   ├── network_parser.jl                    # Network data preparation for SE
│   └── distributions.jl                     # Statistical distribution utilities
└── prob/se.jl                               # State estimation problem definitions
```

**Key functions:**

- `add_measurements!(data, meas)` — attach measurements to network data
- `solve_mc_se(data, form, solver)` — solve multi-conductor state estimation
- `compute_mc_se(data, form, solver)` — compute SE (no optimization)
- `_meas_comp_to_ext_comp(...)` — measurement type conversion between component types
- Chi-squared test (`bad_data/chi_squares_test.jl`) and LNR (`largest_normalized_residuals.jl`) for bad data detection

The `feat/transformer_handling` branch adds active transformer modelling within SE — relevant to TSK-844.

## TSK-844 Active Research Context

**Task**: Timeseries MV/LV transformer tap monitoring via DSSE. Five-phase plan:

1. **Network Setup**: Multinetwork DSS parsing → IVR/EN math model → OLTC configuration → simulate tap change at t=49
2. **Reference Generation**: Timeseries OPF → synthetic measurements with Gaussian noise (τ=0.05) → CSV export
3. **Standard SE**: Run WLS SE → detect anomaly via objective spike at t=49 → localize transformer via residuals
4. **Tap Recovery**: MINLP warm-start SE for discrete tap estimation → validate against ground truth
5. **Visualization**: Pliers-based plots — objective timeseries, residual heatmaps

## Power Systems Domain Capabilities

Technical areas relevant to this codebase:

- Three-phase unbalanced power flow (balanced/unbalanced), sequence components (positive/negative/zero)
- Per-unit system conversions, impedance/admittance matrices, Kron reduction
- Transformer equivalent circuits (T-model, π-model), tap ratio modelling, OLTC/tap changer simulation
- Distribution System State Estimation (DSSE): WLS, WLAV, LAV formulations
- Synthetic measurement generation with configurable noise distributions (Gaussian, etc.)
- Bad data detection: chi-squared hypothesis tests, largest normalized residuals (LNR)
- Measurement placement and observability analysis
- Voltage regulation (OLTC, capacitor banks, DER)

**Key references**: Kersting's *Distribution System Modeling and Analysis*; Abur & Expósito's *Power System State Estimation: Theory and Implementation*

**Standards**: IEC 60076 (power transformers), IEC 61850 (substation communication), IEC 61968/61970 (CIM), IEC 60038 (standard voltages), EN 50160 (voltage quality), IEEE 1547, IEEE C57 series

## PMD Enum Constants — Never Use Bare Names

Pliers does **not** import PMD constants (`MATHEMATICAL`, `ENGINEERING`, `DELTA`, `WYE`, etc.) into its module scope. Using them bare causes `UndefVarError`. Always compare via string:

```julia
# CORRECT
string(data["data_model"]) == "MATHEMATICAL"
string(load["configuration"]) == "DELTA"

# WRONG — UndefVarError at runtime
data["data_model"] == MATHEMATICAL
load["configuration"] == DELTA
```

This applies everywhere in `PMDUtils`, `PMDGraph`, and any other Pliers module.

## Graph Visualization Internals (PMDGraph)

- **Transformer directionality in MATH**: Components named `_virtual_branch..._2` or `_virtual_transformer...2` are reverse-directed. Flip `f_bus`/`t_bus` for names starting with `_virtual_branch`/`_virtual_transformer` AND ending with `_2`/`.2`.
- **`smart_layout()`**: Handles non-tree topologies by temporarily removing transformer edges to compute Buchheim layout, then restores all edges.
- **Duplicate edges**: Use `has_edge()` before adding transformer edges; merge properties rather than adding duplicates (required for `edge_color` and other decorations).
- **Graph type detection**: Use `_is_eng_graph(graph::AbstractMetaGraph)` — do not confuse with `_is_eng(data::Dict)` in PMDUtils.

## Workflow for Significant Tasks

1. **Consult knowledge base**: Read `.github/agents/knowledge/workspace-discovery.md` for prior context and lessons learned.
2. **Plan**: Break task into clear steps. Identify required data, formulations, and code paths.
3. **Research**: Search the workspace, read relevant source files. Use web search for standards/paper references.
4. **Implement**: Write clean, documented Julia code. Follow PMD/PMDSE conventions. Justify assumed values inline.
5. **Validate**: Run code, check results against expected behavior, verify units and magnitudes.
6. **Document**: Record findings, update `.github/agents/knowledge/workspace-discovery.md` with new patterns.

## Engineering Assumptions & Output Standards

When making engineering assumptions, always document them explicitly:

| Parameter | Value | Unit | Source/Justification                         |
|-----------|-------|------|----------------------------------------------|
| (example) | 0.05  | p.u. | Gaussian noise std, typical for smart meters |

- Cite references inline: (Kersting, 2012) or (IEC 60076-1:2011, §4.2)
- Validate result magnitudes (e.g., "voltage at 0.98 p.u. is within EN 50160 ±10% band")

## Constraints

- **DO NOT** silently assume parameter values — always document and justify.
- **DO NOT** skip validation steps even under time pressure.
- **DO NOT** modify PMD or PMDSE source code without explicit user approval — they are upstream packages.
- **DO NOT** use approximate formulations when exact ones are available in PMD/PMDSE.
- **ALWAYS** check which model (ENG vs MATH) you are working with before accessing components.

## Coding Conventions

- Function names: `snake_case`
- **PrettyTables**: Always pass a `DataFrame`, never a matrix with `header` kwarg:
  ```julia
  # CORRECT
  df = DataFrame(Attribute=keys, Value=values)
  pretty_table(df)

  # AVOID
  pretty_table([keys values]; header=["Attribute", "Value"])
  ```
- Docstrings must be compatible with `Documenter.jl`
- Visualization uses the `Makie.jl` ecosystem: `CairoMakie` (static), `GeoMakie` (maps), `GraphMakie` (graphs)
- Examples in `examples/` are converted to documentation tutorials via `Literate.jl`

## Knowledge Base

The authoritative workspace reference and lessons learned are in `.github/agents/knowledge/workspace-discovery.md`. Consult it before starting significant tasks and update it when discovering new patterns or completing significant work.
