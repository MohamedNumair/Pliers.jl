# Pliers.jl Documentation

## Overview

Pliers.jl is a collection of tools that I usually need (like pliers for an electrician) for analyzing power distribution systems. It is designed to be used in conjunction with the [PowerModelsDistribution.jl](https://github.com/lanl-ansi/PowerModelsDistribution.jl) or [PowerModelsDistributionStateEstimation.jl](https://github.com/Electa-Git/PowerModelsDistributionStateEstimation.jl) packages for simplified reporting, analysis and visualization.

## Features

- **Network Visualization**: Plot distribution networks as trees, with coordinates, or on geographic maps
- **Model Exploration**: Generate reports and tables for engineering and mathematical network models
- **State Estimation Tools**: Visualize residuals and analyze measurement data
- **Result Processing**: Transform and analyze power flow results

## Package Organization

Pliers.jl provides all functionality through the main module. For organizational purposes, functions are also re-exported through sub-modules that group related functionality:

### PMDUtils
Re-exports utility functions for PowerModelsDistribution workflows including:
- Solution processing (voltage fluffing, dictification)
- Impedance calculations (Kron reduction, sequence components)
- Network data manipulation
- Engineering and mathematical model exploration

### PMDSEUtils
Re-exports utility functions for PowerModelsDistributionStateEstimation workflows including:
- State estimation result visualization
- Measurement residual analysis
- Measurement data processing and writing

### PMDGraph
Re-exports plotting functions for network visualization including:
- Network tree visualization
- Coordinate-based network plotting
- Geographic map overlays
- Bus phasor diagram plotting

## Installation

Currently Pliers is still under development and not registered. You can install it by cloning the repository and adding it to your Julia environment:

```julia
using Pkg
Pkg.develop(url="https://github.com/MohamedNumair/Pliers.jl")
```

## Quick Start

```julia
using Pliers

# Set journal-quality theme for plots
inch, pt, cm, ieeecolumn, ieee2column = set_journal_theme()

# Generate a report for an engineering model
# eng_report(eng)

# Visualize network as a tree
# plot_network_tree(network_graph)
```

## Author

This package was developed by Mohamed Numair ([mnumair.com](https://mnumair.com))

