module GravityModels

using LinearAlgebra

abstract type AbstractGravity end

export GravityModel, gravity, gravity!
include("keplerian.jl")
include("harmonics.jl")

export geodetic2geocentric, geodetic2geocentric!, gravityMap
include("utils.jl")

end
