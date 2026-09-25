
using SatelliteToolboxGravityModels
const STB = SatelliteToolboxGravityModels.GravityModels

using GravityModels
using BenchmarkTools

egm2008 = STB.load(IcgemFile, fetch_icgem_file(:EGM2008))
workspace = STB.Workspace(egm2008, max_degree=100)
gm = GravityModel(GravityModels.EGM2008)
g = zeros(3)

err = -1.0
for i in 1:1000
    pos = normalize(randn(3)) * 7000.0e3 + randn(3) .* 100e3
    gtest = STB.gravitational_acceleration(egm2008, pos; workspace, max_degree=100)
    gravity!(gm, pos, g)
    norm(g - gtest)
    err = max(err, maximum(abs, g - gtest))
end
@show err

# @btime STB.gravitational_acceleration(egm2008, $pos; workspace, max_degree=100)
# @btime gravity!($gm, $pos, $g)
