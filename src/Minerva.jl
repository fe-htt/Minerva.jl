"""
    FuelCycleSimulation

A package for simulating fusion reactor fuel cycle using ModelingToolkit.
Provides simulation orchestration, ensemble analysis, and postprocessing tools
for tritium fuel cycle models.
"""
module Minerva

# Re-export key dependencies users will need
# using Reexport

# Core simulation infrastructure
#@reexport using DifferentialEquations
#@reexport using ModelingToolkit
#@reexport using ModelingToolkitStandardLibrary
#@reexport using FuelCycleLibrary
#@reexport using SciCompDSL

# Essential for user input files
#@reexport using FuelCycleLibrary.PlasmaFunctions
#@reexport using FuelCycleLibrary.Base.Atomic
#@reexport using ModelingToolkit: t_nounits as t, D_nounits as D
#@reexport using ModelingToolkitStandardLibrary.Blocks




using Printf
using IterTools
using CSV
using DataFrames

# Postprocessing



include("simulation_engine.jl")
#include("postprocessing/drawio_utils.jl")
#include("postprocessing/inventory_viz.jl")
#include("postprocessing/plotting_utils.jl")



# Public API
export process_sim_settings, run_ensemble_sim, setup_ensemble_problem,
    save_sim_results, print_sim_header

end
