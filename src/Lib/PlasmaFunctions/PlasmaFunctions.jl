"""
Library of basic models.
"""
module PlasmaFunctions
using QuadGK


export get_nProfile, get_TProfile, get_nAveraged, get_nAveraged_radius, get_nAveraged_disk, 
       get_nAveraged_radius_2018, get_nAveraged_disk_2018, get_nAveraged2

export get_nAvg_powerCorrectionFactor, get_nAvg_radius_powerCorrectionFactor, get_nAvg_disk_powerCorrectionFactor, 
       get_nAvg_radius_powerCorrectionFactor_2018, get_nAvg_disk_powerCorrectionFactor_2018

export get_reactionRates, get_DT_reactionRate

export saveSim



include("functions.jl")

end