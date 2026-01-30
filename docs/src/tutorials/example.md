```@meta
EditURL = "../../../examples/example.jl"
```

````@example example
using Pkg
Pkg.add(PackageSpec(path=pwd()))

using Pliers
using PowerModelsDistribution
````

parse dss file

````@example example
dss_file = joinpath(@__DIR__, "..", "examples", "IEEE_13_Node_Test_Feeder.dss")
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

