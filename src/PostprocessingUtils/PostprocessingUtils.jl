"""
Utility functions for postprocessing simulation results.
"""
module PostprocessingUtils
using CSV, DataFrames, LightXML, Colors, Printf

export update_drawio_labels_with_modifications2
export set_drawio_labels_to_na
export update_drawio_inventory_colours

include("drawio_utils.jl")
include("plotting_utils.jl")
include("inventory_viz.jl")

end #module