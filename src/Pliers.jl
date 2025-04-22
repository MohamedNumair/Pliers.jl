module Pliers


# data structure
using DataFrames
using CSV
using FileIO
using Serialization
# pretty terminal packages
using Crayons
using Crayons.Box
using PrettyTables

# plotting packages
using Makie
using CairoMakie
using WGLMakie
if Sys.iswindows()
    using GLMakie
end
using GraphMakie
using GeoMakie
using Proj
using Tyler
using Tyler.TileProviders
using Tyler.MapTiles
using Tyler.Extents

# data analysis packages
using Statistics
using LinearAlgebra
using Graphs
using MetaGraphs
# using Dates
# using StatsPlots
# using StatsBase
# using Distributions
# using Random  

# Power Distribution Tools
using PowerModelsDistribution
using PowerModelsDistributionStateEstimation
using Ipopt

# pkg const
const pkg_name = "Pliers"
const BASE_DIR = dirname(@__DIR__)

_N_IDX = 4 # representing the neutral index in terminals 

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

author() = println("This package was developped by Mohamed Numair ✪ ω ✪")


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

#PMD
include("core/PMD/pmd_utils.jl")
include("core/PMD/network_graph.jl")
include("core/PMD/results_explorer.jl")
include("core/PMD/eng_explorer.jl")
include("core/PMD/math_explorer.jl")
include("core/PMD/network_plotting.jl")

#PMDSE
include("core/PMDSE/pmdse_utils.jl")


end # module Pliers
