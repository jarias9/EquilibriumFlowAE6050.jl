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
include("reactionmodels.jl")
include("reactionsystem.jl")
include("equilibriumconstant.jl")
include("calculateEquilibrium.jl")
include("constants.jl")
include("partitionfunctions.jl")
include("equilibriumTP.jl")
include("equilibriumHP.jl")
include("equilibriumentropy.jl")
include("flowmodels.jl")
include("equilibriumflow.jl")

# Export relevant methods
export equilibriumTP, equilibriumHP, equilibriumflow

end
