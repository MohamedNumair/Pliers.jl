# Pliers.jl Copilot Instructions

## Project Context
Pliers.jl is a Julia package providing tools for analyzing power distribution systems, designed to work with `PowerModelsDistribution.jl` (PMD) and `PowerModelsDistributionStateEstimation.jl` (PMDSE). It focuses on reporting, analysis, and visualization.

## Architecture
- **Core Logic**: Located in `src/`.
- **Submodules**: Functionality is organized into submodules under `src/modules/` (e.g., `PMDUtils`, `PMDGraph`) and re-exported by `src/Pliers.jl`.
- **Data Models**: Engineering data and mathematical models are represented as `Dict{String, Any}`.
- **Visualization**: Uses `Makie.jl` ecosystem (`CairoMakie`, `GeoMakie`, `GraphMakie`).

## Coding Conventions
- **DataFrames & Tables**: 
    - When using `PrettyTables.jl`, **always prefer passing a `DataFrame`** with defined column names rather than passing a matrix and using the `header` keyword argument. The `header` keyword is not reliably supported for all matrix inputs in the current configuration.
    - **CORRECT**: 
      ```julia
      df = DataFrame(Attribute=keys, Value=values)
      pretty_table(df)
      ```
    - **AVOID**: 
      ```julia
      pretty_table([keys values]; header=["Attribute", "Value"])
      ```
- **Functions**: Use `snake_case` for function names.
- **Documentation**: Use logical docstrings compatible with `Documenter.jl`.

## Development Workflows
- **Environment**: Always activate the package environment before running code: `julia --project=.` or `pkg> activate .`.
- **Testing**: Run tests via `julia --project=. test/runtests.jl` or `pkg> test`.
- **Documentation**: Build docs using `julia --project=docs docs/make.jl`. The documentation build uses `Literate.jl` to generate tutorials from `examples/`.

## Key Dependencies
- `PowerModelsDistribution`: Implicit dependency for data structures.
- `DataFrames`: For tabular data handling.
- `CairoMakie`: For static plotting.
- `PrettyTables`: For terminal output of tables.
