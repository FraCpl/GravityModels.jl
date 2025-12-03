struct GravityKeplerian{T} <: AbstractGravity
    μ::T
end

Base.broadcastable(x::GravityKeplerian) = Ref(x)

function GravityModel(μ::T) where {T}
    return GravityKeplerian(μ)
end

@inline function gravity!(μ::T, pos::AbstractVector{T}, g::AbstractVector{T}) where {T}
    x, y, z = pos
    r2 = x*x + y*y + z*z
    r3 = r2*sqrt(r2)
    c = -μ/r3
    g[1] = c*x;
    g[2] = c*y;
    g[3] = c*z
    return nothing
end

@inline function gravity!(GH::GravityKeplerian, pos::AbstractVector{T}, g::AbstractVector{T}) where {T}
    gravity!(GH.μ, pos, g)
    return nothing
end

@inline function gravity(GH::GravityKeplerian, pos::AbstractVector{T}) where {T}
    g = zero(pos)
    gravity!(GH.μ, pos, g)
    return g
end
