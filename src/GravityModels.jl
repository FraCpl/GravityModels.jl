module GravityModels

using LinearAlgebra

abstract type AbstractGravity end

const DATA::String = replace(pathof(GravityModels)[1:(end - 20)], "\\" => "/")*"data/"
const EGM2008::String = DATA*"earth/EGM2008.gfc"
const JMM3::String = DATA*"mars/gmm3_120.sha.gfc"
const JGMRO::String = DATA*"mars/jgmro_120f_sha.gfc"
const GL0660B::String = DATA*"moon/GL0660B.gfc"

export GravityModel, gravity, gravity!
include("keplerian.jl")
include("harmonics.jl")

export geodetic2geocentric, geodetic2geocentric!, gravityMap
include("utils.jl")

end
