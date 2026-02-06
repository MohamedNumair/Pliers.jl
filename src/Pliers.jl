"""
    Pliers

A Julia package providing tools for analyzing power distribution systems. Designed to be used
in conjunction with PowerModelsDistribution.jl and PowerModelsDistributionStateEstimation.jl
for simplified reporting, analysis, and visualization.

# Sub-modules
The package provides optional sub-module access for organized imports:
- `PMDUtils`: Re-exports PMD-related utility functions
- `PMDSEUtils`: Re-exports PMDSE-related utility functions  
- `PMDGraph`: Re-exports plotting functions

# Author
Mohamed Numair (mnumair.com)
"""
module Pliers

# data structure
using DataFrames
using CSV
using FileIO
using Serialization
using Suppressor

# pretty terminal packages
using Crayons
using Crayons.Box
using PrettyTables

# plotting packages
using Makie
# using MakieCore
using CairoMakie
using WGLMakie
if Sys.iswindows()
    # using GLMakie
end

# data analysis packages
using Statistics
using LinearAlgebra


# pkg const
const pkg_name = "Pliers"
const BASE_DIR = dirname(@__DIR__)


# data structure
const _DF = DataFrames
const _CSV = CSV

# pretty terminal packages
const _CRN = Crayons
const _PT = PrettyTables


# plotting packages
const _MK = Makie
const _CMK = CairoMakie
const _WGLMK = WGLMakie
# const _GLMK = GLMakie

# https://github.com/PumasAI/SummaryTables.jl?  # SummaryTable.jl is amazing in outputting tables right away

# data analysis packages
const _STAT = Statistics

# Power Distribution Tools
# const _PMD = PowerModelsDistribution
# const _PMDSE = PowerModelsDistributionStateEstimation

"""
    author()

Print information about the package author.
"""
author() = println("This package was developed by Mohamed Numair (mnumair.com)")


# included functionalities


# helper functions
include("core/styles.jl")
include("core/utils.jl")

#IO
include("core/export.jl")

include("io/pliers-io.jl")
include("io/networks_io.jl")
include("io/util-enwl-networks.jl")
include("io/util-spanish-networks.jl")


# Sub-modules for organized namespace
include("modules/PMDUtils.jl")
include("modules/PMDSEUtils.jl")
include("modules/PMDGraph.jl")

end # module Pliers
