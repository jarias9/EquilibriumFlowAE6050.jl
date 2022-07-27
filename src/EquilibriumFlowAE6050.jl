module EquilibriumFlowAE6050

# Declare Dependencies
using LinearAlgebra
using Plots
using Printf
using Unitful
using Parameters
using UnitfulRecipes
using ForwardDiff
using LeastSquaresOptim
using Roots

# Import scripts
include("flowmodels.jl")
include("reactionmodels.jl")
include("reactionsystem.jl")
include("equilibriumconstant.jl")
include("calculateEquilibrium.jl")
include("constants.jl")
include("partitionfunctions.jl")
include("equilibriumTP.jl")
include("equilibriumHP.jl")
include("equilibriumentropy.jl")
include("equilibriumflow.jl")

# Export constants 
export constants

# Export models and flow problems
export Air5s, Air7s, Air11s, Air13s, CO2_6s
export NormalShock, ObliqueShock, DeflectionAngle, ShockAngle

# Export relevant methods
export equilibriumTP, equilibriumHP, equilibriumflow

end
