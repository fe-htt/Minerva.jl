module Lib

Na::Float64 = 6.02214076e+23 # Avogadro number [mol-1]
Rg::Float64 = 8.314 # Ideal gas costant [J mol-1 K-1]
λT::Float64 = 1.782607697460566e-9 # tritium decay rate [s-1]

include("Base/Base.jl")
include("PlasmaFunctions/PlasmaFunctions.jl")


export Na, Rg, λT
end

