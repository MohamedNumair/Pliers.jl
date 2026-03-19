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

# PMD Data Models (ENG vs MATH)
When translating code between Engineering (ENG) and Mathematical (MATH) models, observe the following structural differences:

- **Component Mapping**:
    - **Lines**: Represented as `line` in ENG, but mapped to `branch` in MATH.
    - **Transformers**: Represented as `transformer` in ENG. In MATH, a transformer is decomposed into **two ideal transformers and three virtual branches**.
- **Topological Artifacts**:
    - **Virtual Elements**: Explicit `_virtual_bus` and `_virtual_branch` elements are often added in the MATH model (e.g., connecting voltage sources or decomposing transformers) to satisfy mathematical formulation requirements.
    - **Indexing**: ENG models use arbitrary String IDs. MATH models use Integer indexing.

## Graph Visualization & PMD Internals (`PMDGraph`)
- **Transformer Decomposition in MATH Models**:
    - Transformers in MATH models are decomposed into ideal transformers and virtual branches.
    - **Directionality Issue**: Components like `_virtual_branch..._2` and `_virtual_transformer...2` are often defined in the reverse direction of the logical "flow" needed for tree plotting.
    - **Fix**: When constructing the graph, check for names starting with `_virtual_branch`/`_virtual_transformer` AND ending with `_2`/`.2`. Flip the `f_bus` and `t_bus` for these edges to maintain a consistent rooted tree structure.
- **Layout Algorithms**:
    - **`smart_layout`**: Robustly handles non-tree graphs by temporarily removing "logical" edges (like transformers that create cycles) to compute a `Buchheim` tree layout, then applies positions to the full graph.
- **Edge Properties**:
    - When adding transformers to the graph, use `has_edge` to check for existing branches (e.g., lines) between the same nodes. **Merge properties** instead of adding duplicate edges to ensure plotting decorations (like `edge_color`) work correctly.
    - Use `_is_eng_graph(graph)` for `AbstractMetaGraph` inputs to avoid MethodError conflicts with `PMDUtils._is_eng(data::Dict)`.