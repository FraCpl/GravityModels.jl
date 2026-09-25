struct GravityHarmonics <: AbstractGravity
    # Model parameters
    μ::Float64
    order::Int
    Rref::Float64
    SCALING::Float64
    S::Matrix{Float64}
    C::Matrix{Float64}

    # Precomputed variables
    gnmOj::Vector{Float64}
    hnmOj::Vector{Float64}
    enm::Vector{Float64}
    sectorial::Vector{Float64}

    # Allocations
    aOrN::Vector{Float64}
    cosλ::Vector{Float64}
    sinλ::Vector{Float64}
    pnm0Plus2::Vector{Float64}
    pnm0Plus1::Vector{Float64}
    pnm0::Vector{Float64}
    pnm1::Vector{Float64}
end

Base.broadcastable(x::GravityHarmonics) = Ref(x)

function GravityModel(coeffFile::String; order::Int64=-1, SCALING=0.0)
    # Only compatible with STATIC gravity field models downloaded from
    # http://icgem.gfz-potsdam.de/home

    # Parse model metadata
    μ = 0.0
    R = 0.0
    maxOrder = 0
    open(coeffFile, "r") do file
        for ln in eachline(file)
            data = split(ln)
            isempty(data) && continue

            if occursin("gravity_constant", data[1])
                μ = parse(Float64, data[2])

            elseif occursin("radius", data[1])
                R = parse(Float64, data[2])

            elseif occursin("max_degree", data[1])
                maxOrder = parse(Int, data[2])
                break
            end
        end
    end
    order = order < 0 ? maxOrder : min(order, maxOrder)
    @assert order > 0

    # Read coefficients
    C = zeros(Float64, order + 1, order + 1)
    S = zeros(Float64, order + 1, order + 1)
    open(coeffFile, "r") do file
        for ln in eachline(file)
            data = split(ln)
            isempty(data) && continue
            data[1] != "gfc" && continue

            n = parse(Int, data[2])
            n > order && continue

            m = parse(Int, data[3])
            m > order && continue

            C[n + 1, m + 1] = parse(Float64, data[4])
            S[n + 1, m + 1] = parse(Float64, data[5])
            n == order && m == order && break
        end
    end
    C[1, 1] = 1.0

    # ================= PRECOMPUTATIONS NEW ALGORITHM ======================= #
    # pre-compute the recursion coefficients corresponding to equations 19 and 22
    # from Holmes and Featherstone paper
    # for cache efficiency, elements are stored in the same order they will be used
    # later on, i.e. from rightmost column to leftmost column
    degree = order
    N = degree * (degree + 1) ÷ 2 - 1
    gnmOj = Vector{Float64}(undef, N)
    hnmOj = Vector{Float64}(undef, N)
    enm = Vector{Float64}(undef, N)
    i = 1
    for m in degree:-1:0
        j = 2.0 - 1.0*(m > 0)
        for n in max(2, m + 1):degree
            f = (n - m)*(n + m + 1)
            gnmOj[i] = 2*(m + 1) / sqrt(j*f)
            hnmOj[i] = sqrt((n + m + 2)*(n - m - 1) / (j*f))
            enm[i] = sqrt(f / j)
            i += 1
        end
    end

    # scaled sectorial terms corresponding to equation 28 in Holmes and Featherstone paper
    sectorial = zeros(degree + 1)
    sectorial[1] = scalb(1.0, -SCALING)
    sectorial[2] = sqrt(3.0)*sectorial[1]
    for m in 2:degree
        sectorial[m + 1] = sqrt((2*m + 1) / (2*m))*sectorial[m]
    end

    cosλ = zeros(degree+1); cosλ[1] = 1.0
    return GravityHarmonics(μ, order, R, SCALING, S, C, gnmOj, hnmOj, enm, sectorial,
        ones(degree+1), cosλ, zeros(degree+1), zeros(degree+1), zeros(degree+1),
        zeros(degree+1), zeros(degree+1))
end

Base.copy(g::GravityHarmonics) = GravityHarmonics(
    g.μ,
    g.order,
    g.Rref,
    g.SCALING,
    copy(g.S),
    copy(g.C),
    copy(g.gnmOj),
    copy(g.hnmOj),
    copy(g.enm),
    copy(g.sectorial),
    copy(g.aOrN),
    copy(g.cosλ),
    copy(g.sinλ),
    copy(g.pnm0Plus2),
    copy(g.pnm0Plus1),
    copy(g.pnm0),
    copy(g.pnm1),
)

@inline scalb(x, n) = x * exp2(n)

#  gravity(GH, x, y, z)
#  This function computes the gravitational acceleration of a gravitational
#  potential developed in spherical harmonics. This function is much faster
#  than the equivalent one in STAR, and robust to very high order/degree.
#
#  INPUTS:
#  GH                   Gravity Harmonics data structure (to be obtained with GravityModel(...))
#  x, y, z  [m]         Position of the vehicle in planet centered, planet fixed reference frame
#
#  OUTPUTS:
#  gx, gy, gz  [m/s^2]  Gravity acceleration vector in planet fixed reference frame
#
# Reference:
# [1] S. A. Holmes, W. E. Featherstone, A unified approach to the Clenshaw
# summation and the recursive computation of very high degree and order
# normalised associated Legendre functions. Journal of Geodesy (2002) 76:
# 279–299. DOI 10.1007/s00190-002-0216-2.
#
# Adapted from OREKIT
function gravity(GH::GravityHarmonics, x::T, y::T, z::T) where {T}
    # Extract data from GH
    pnm0      = GH.pnm0
    pnm1      = GH.pnm1
    pnm0Plus1 = GH.pnm0Plus1
    pnm0Plus2 = GH.pnm0Plus2
    order     = GH.order
    Rref      = GH.Rref
    μ         = GH.μ
    S         = GH.S
    C         = GH.C
    aOrN      = GH.aOrN
    cosλ      = GH.cosλ
    sinλ      = GH.sinλ
    SCALING   = GH.SCALING

    # Compute polar coordinates
    ρ2 = x * x + y * y
    r2 = ρ2 + z * z
    r = sqrt(r2)
    ρ = sqrt(ρ2)
    t = z / r         # cos(theta), where theta is the polar angle
    u = ρ / r         # sin(theta), where theta is the polar angle
    tOu = z / max(ρ, eps(T))

    createDistancePowersArray!(GH, Rref / r)    # compute distance powers
    createCosSinArrays!(GH, x / ρ , y / ρ)      # compute longitude cosines/sines

    # outer summation over order
    index = 1
    U = zero(T)
    dUx = zero(T)
    dUy = zero(T)
    dUz = zero(T)
    # fill!(GH.pnm0, 0.0)
    # fill!(GH.pnm1, 0.0)
    fill!(pnm0Plus1, 0.0)
    fill!(pnm0Plus2, 0.0)

    @inbounds for m in order:-1:0
        # compute tesseral terms with derivatives
        index = computeTesseral!(GH, m, index, t, u, tOu, pnm0, pnm1, pnm0Plus1, pnm0Plus2)

        # compute contribution of current order to field (equation 5 of the paper)
        # inner summation over degree, for fixed order
        sumDegreeS = 0.0
        sumDegreeC = 0.0
        dSumDegreeSdR = 0.0
        dSumDegreeCdR = 0.0
        dSumDegreeSdTheta = 0.0
        dSumDegreeCdTheta = 0.0
        mp1 = m + 1

        @inbounds for n in max(2, m):order
            np1 = n + 1
            qSnm = aOrN[np1] * S[np1, mp1]
            qCnm = aOrN[np1] * C[np1, mp1]
            nOr = n / r
            s0 = pnm0[np1] * qSnm
            c0 = pnm0[np1] * qCnm
            sumDegreeS += s0
            sumDegreeC += c0
            dSumDegreeSdR -= nOr * s0
            dSumDegreeCdR -= nOr * c0
            dSumDegreeSdTheta += pnm1[np1] * qSnm
            dSumDegreeCdTheta += pnm1[np1] * qCnm
        end

        # Contribution to outer summation over order
        smλ = sinλ[mp1]
        cmλ = cosλ[mp1]
        U = muladd(U, u, muladd(smλ, sumDegreeS, cmλ * sumDegreeC))
        dUx = dUx * u + smλ * dSumDegreeSdR + cmλ * dSumDegreeCdR
        dUy = dUy * u + smλ * dSumDegreeSdTheta + cmλ * dSumDegreeCdTheta
        dUz = dUz * u + m * (cmλ * sumDegreeS - smλ * sumDegreeC)

        # Rotate the recursion arrays
        pnm0Plus2, pnm0Plus1, pnm0 = pnm0Plus1, pnm0, pnm0Plus2
    end

    # Scale back
    if SCALING > 0
        U = scalb(U, SCALING)
        dUx = scalb(dUx, SCALING)
        dUy = scalb(dUy, SCALING)
        dUz = scalb(dUz, SCALING)
    end

    # apply the global mu/r factor
    muOr = μ / r
    dUx = dUx * muOr - U * muOr / r
    dUy *= -muOr
    dUz *= muOr

    # Convert Gradient from Spherical to Cartesian Coordinates and add C[1, 1] = C₀₀ term
    rI = x / r
    rJ = y / r
    c00 = -μ / (r2 * r)# * C[1, 1]      # By construction C[1, 1] = 1
    gx = rI * dUx - rI * t / ρ * dUy - rJ * r / ρ2 * dUz + c00 * x
    gy = rJ * dUx - rJ * t / ρ * dUy + rI * r / ρ2 * dUz + c00 * y
    gz = t * dUx + ρ / r2 * dUy + c00 * z

    return gx, gy, gz
end

# This function computes normalized associated legendre functions
function computeTesseral!(GH::GravityHarmonics, m, index, t, u, tOu, pnm0, pnm1, pnm0Plus1, pnm0Plus2)
    order = GH.order
    gnmOj = GH.gnmOj
    hnmOj = GH.hnmOj
    enm = GH.enm
    nmax = max(2, m)

    # initialize recursion from sectorial terms
    n = nmax
    if n == m
        @inbounds pnm0[n + 1] = GH.sectorial[n + 1]
        n += 1
    end

    # compute tesseral values
    localIndex = index
    u2 = u * u
    mtOu = m * tOu

    @inbounds for k in (n + 1):(order + 1)
        # value (equation 27 of the paper)
        pnm0[k] = gnmOj[localIndex] * t * pnm0Plus1[k] - hnmOj[localIndex] * u2 * pnm0Plus2[k]
        localIndex += 1
    end

    # initialize recursion from sectorial terms
    n = nmax
    if n == m
        @inbounds pnm1[n + 1] = mtOu * pnm0[n + 1]
        n += 1
    end

    # compute tesseral values and derivatives with respect to polar angle
    localIndex = index
    @inbounds for k in (n + 1):(order + 1)
        # first derivative (equation 30 of the paper)
        pnm1[k] = mtOu * pnm0[k] - enm[localIndex] * u * pnm0Plus1[k]
        localIndex += 1
    end
    return localIndex
end

# This function computes (Rref/r)^n
@inline function createDistancePowersArray!(GH::GravityHarmonics, aOr)
    GH.aOrN[2] = aOr
    @inbounds for n in 2:GH.order
        p = fld(n, 2) + 1
        q = n - p + 2
        GH.aOrN[n + 1] = GH.aOrN[p] * GH.aOrN[q]
    end
    return nothing
end

# This function computes cos(m*lon) and sin(m*lon)
@inline function createCosSinArrays!(GH::GravityHarmonics, cosLambda, sinLambda)
    GH.cosλ[2] = cosLambda
    GH.sinλ[2] = sinLambda

    # fill up array
    @inbounds for m in 2:GH.order
        # m*lambda is split as p*lambda + q*lambda, trying to avoid
        # p or q being much larger than the other. This reduces the number of
        # intermediate results reused to compute each value, and hence should limit
        # as much as possible roundoff error accumulation
        # (this does not change the number of floating point operations)
        p = fld(m, 2) + 1
        q = m - p + 2

        GH.cosλ[m + 1] = GH.cosλ[p] * GH.cosλ[q] - GH.sinλ[p] * GH.sinλ[q]
        GH.sinλ[m + 1] = GH.sinλ[p] * GH.cosλ[q] + GH.cosλ[p] * GH.sinλ[q]
    end

    return nothing
end
