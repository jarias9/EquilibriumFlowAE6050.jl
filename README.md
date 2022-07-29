# EquilibriumFlowAE6050

Library implementing an Equilibrium Flow solver as part of the Final Project for AE 6050 for Summer 2022.

# Simple Example

The library has three main functions, with two functions used for calculating the equilibrium state of a mixture of gases, and another function for solving equilibrium flow problems. The syntax for each is as follows 

First initialize input parameters. Note: the solvers expect inputs with units via Unitful.jl
```julia
T = 3000u"K" # 3000 K
P = 101325u"Pa" # 101325 Pa, 1 atm
```

The equilibrium composition solvers can be accessed via
```julia
using EquilibriumFlowAE6050
# T, P solver
h_mix, rho_mix, s_mix, X = equilibriumTP(Air11s(), T, P)

# h, P solver
T, rho_mix, s_mix, X = equilibriumHP(model::GasModel, h_mix, P)
```
Four air models are included being ```Air5s()```, ```Air7s()```, ```Air11s()```, ```Air13s()```.

To solve a flow problem the ```equilibriumflow(problem, model, ...)``` function is used. There are two types of problems that can currently be solved, normal shocks, and oblique shocks. An example of the interface for each is as follows
```julia
# Solve Normal Shock Problem
T2, h2, P2, rho2, u2, s2, X2 = equilibriumflow(NormalShock(), Air11s(), u1, T1, P1)
# Solve Oblique Shock Problem, Known Shock Wave Angle
T2, h2, P2, rho2, u2, s2, X2, theta = equilibriumflow(ObliqueShock(), ShockAngle(), Air11s(), angle, u1, T1, P1)
# Solve Oblique Shock Problem, Known Deflection Angle
T2, h2, P2, rho2, u2, s2, X2, beta = equilibriumflow(ObliqueShock(), DeflectionAngle(), Air11s(), angle, u1, T1, P1)
```
