using GravityModels
using Test
using LinearAlgebra
using DelimitedFiles

function TEST_gravity()
    # Ground truth data come from: https://icgem.gfz-potsdam.de/calcgrid
    # Functional selection: gravitation_ell, WGS84 reference system, 100x100 order/degree
    data = readdlm(
        dirname(pathof(GravityModels))[1:(end-3)]*"test\\gravity_earth_EGM2008_100x100.gdf",
        skipstart = 34,
    )
    GM = GravityModel(GravityModels.EGM2008)

    errPerc = -1.0;
    pos = zeros(3);
    grav = zeros(3)
    for r in eachrow(data)
        lond, latd, g = r
        geodetic2geocentric!(pos, latd*π/180, lond*π/180, 0.0)
        gravity!(GM, pos, grav)
        errPerc = max(errPerc, abs(1.0 - norm(grav)/(1e-5g)))
    end

    return errPerc
end

@testset "GravityModels.jl" begin
    @test TEST_gravity() < 1e-9
end
