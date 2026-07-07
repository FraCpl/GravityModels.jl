struct GravityKeplerian{T} <: AbstractGravity
    μ::T
end

Base.broadcastable(x::GravityKeplerian) = Ref(x)

GravityModel(μ::T) where {T} = GravityKeplerian(μ)

@inline gravity(GH::GravityKeplerian{M}, x::T, y::T, z::T) where {M, T} = gravity(GH.μ, x, y, z)

@inline function gravity(μ::M, x::T, y::T, z::T) where {M, T}
    r2 = x*x + y*y + z*z
    rinv = 1 / sqrt(r2)
    rinv3 = rinv * rinv * rinv
    c = -μ * rinv3
    return c * x, c * y, c * z
end
