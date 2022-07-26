using EquilibriumFlowAE6050
using Unitful

T1 = 3000u"K"
P1 = 101325u"Pa"

h, rho, s, X = equilibriumTP(Air11s(), T1, P1)