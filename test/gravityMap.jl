using GravityModels
using LinearAlgebra
using GLMakie

GM = GravityModel(GravityModels.GL0660B, order=60)
# GM = GravityModel(GravityModels.EGM2008, order=60)

# GM.C[1, 1] = 0.0    # No μ
# GM.C[3, 1] = 0.0    # No J2

function plotGravityMap(GH::GravityModels.GravityHarmonics, altitude=100e3; levels=100, N=400)
    lat, lon, g = gravityMap(GH, altitude; N=N)

    g0 = GH.μ/(GH.Rref + altitude)^2       # Keplerian gravity at the same input altitude
    f = Figure(); display(f)
    ax1 = Axis(f[1, 1]; aspect=DataAspect(), xlabel="Longitude [deg]", ylabel="Latitude [deg]", title="Gravity magnitude at $(altitude/1000) km")
    co = contourf!(ax1, lon, lat, (norm.(g) .- g0)./g0*100; levels=levels, colormap=:grays1)
    Colorbar(f[1, 2], co, label = "Δg [% of μ/r²]")
end

plotGravityMap(GM)
