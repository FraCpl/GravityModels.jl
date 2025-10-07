using GravityModels
using BenchmarkTools
using LinearAlgebra

function gr()
    GM = GravityModel(dirname(pathof(GravityModels))[1:end-3]*"data\\earth\\EGM2008.gfc")
    p = [6700e3*normalize(randn(3)) for _ in 1:10000]
    g = zero(p);

    # @profview gravity!.(GM, p, g)
    @btime gravity!.($GM, $p, $g)
    return
end

gr()
