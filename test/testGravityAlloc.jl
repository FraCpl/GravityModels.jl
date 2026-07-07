using GravityModels
using BenchmarkTools
using LinearAlgebra

function gr()
    GM = GravityModel(dirname(pathof(GravityModels))[1:(end - 3)]*"data\\earth\\EGM2008.gfc")
    x, y, z = 6700e3*normalize(randn(3))
    p = [x; y; z]
    g = zero(p)

    @btime gravity($GM, $x, $y, $z)
    @btime gravity!($GM, $p, $g)
    return nothing
end

gr()
