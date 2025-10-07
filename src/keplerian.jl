struct GravityKeplerian{T} <: AbstractGravity
    μ::T
end

Base.broadcastable(x::GravityKeplerian) = Ref(x)

function GravityModel(μ::Real)
    return GravityKeplerian(μ)
end

@inline function gravity!(μ::T, pos::Vector{T}, g::Vector{T}) where T
    x, y, z = pos
    r = sqrt(x*x + y*y + z*z)
    c = -μ/(r*r*r)
    g[1] = c*x; g[2] = c*y; g[3] = c*z
    return
end

@inline function gravity!(GH::GravityKeplerian, pos::Vector{T}, g::Vector{T}) where T
    gravity!(GH.μ, pos, g)
    return
end

@inline function gravity(GH::GravityKeplerian, pos::Vector{T}) where T
    g = Vector{T}(undef, 3)
    gravity!(GH.μ, pos, g)
    return g
end
