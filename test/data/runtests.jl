# Loading Packages
using PowerModelsDistribution
using Crayons.Box
using DataFrames
using PrettyTables
using Graphs
using GraphMakie
using MetaGraphs
using WGLMakie
using GLMakie
using FileIO
using CSV
using CairoMakie
using Tyler
using Proj


using Tyler.TileProviders
using Tyler.MapTiles
using Tyler.Extents

using Pliers
using Test
# Loading test data
pmd_path = joinpath(dirname(pathof(PowerModelsDistribution)), "..")
pmd_test_models = [
    joinpath(pmd_path, "test/data/opendss/case3_balanced.dss") => "case3_balanced",
    joinpath(pmd_path, "test/data/opendss/case3_unbalanced.dss") => "case3_unbalanced",
    joinpath(pmd_path, "test/data/opendss/case3_balanced_battery.dss") => "case3_balanced_battery",
    joinpath(pmd_path, "test/data/opendss/case5_phase_drop.dss") => "case5_phase_drop",
    joinpath(pmd_path, "test/data/opendss/ut_trans_2w_yy_oltc.dss") => "ut_trans_2w_yy_oltc",
    joinpath(pmd_path, "test/data/opendss/case3_balanced_battery.dss") => "case3_balanced_battery",
]
en_cases_files = readdir(joinpath(pmd_path, "test/data/en_validation_case_data"); join=true)


