using GravityModels
using BenchmarkTools
using LinearAlgebra
using ForwardDiff

@warn "not working with Gravity Harmonics"
GM = GravityModel(dirname(pathof(GravityModels))[1:(end - 3)]*"data\\earth\\EGM2008.gfc")
# GM = GravityModel(3.987e14)

p = 6700e3*normalize(randn(3))
g = zero(p)

f!(y, x) = gravity!(GM, x, y)
J = zeros(3, 3)
ForwardDiff.jacobian!(J, f!, g, p)


f!(g, p)
