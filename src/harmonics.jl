struct GravityHarmonics <: AbstractGravity
    μ::Float64
    order::Int
    Rref::Float64
    S::Matrix{Float64}
    C::Matrix{Float64}
    SCALING::Float64
    gnmOj::Vector{Float64}
    hnmOj::Vector{Float64}
    enm::Vector{Float64}
    sectorial::Vector{Float64}
    dU::Vector{Float64}
    aOrN::Vector{Float64}
    cosλ::Vector{Float64}
    sinλ::Vector{Float64}
    T::Matrix{Float64}
    pnm0Plus2::Vector{Float64}
    pnm0Plus1::Vector{Float64}
    pnm0::Vector{Float64}
    pnm1::Vector{Float64}
end

Base.broadcastable(x::GravityHarmonics) = Ref(x)

function GravityModel(coeffFile::String; order::Int64=-1)

    SCALING = 0.0

    # Only compatible with STATIC gravity field models downloaded from
    # http://icgem.gfz-potsdam.de/home
    μ = 0.0; R = 0.0
    C = zeros(1, 1); S = zeros(1, 1)
    file = open(coeffFile, "r")
    for ln in eachline(file)
        data = split(ln," "; keepempty=false)
        if length(data) == 1
            data = split(ln,"\t"; keepempty=false)
        end
        if length(data) > 1
            if occursin("gravity_constant", data[1])
                μ = parse(Float64, data[2])
            elseif data[1] == "radius"
                R = parse(Float64, data[2])
            elseif data[1] == "max_degree"
                maxOrder = parse(Int64, data[2])
                if order < 0
                    order = maxOrder
                else
                    order = min(order, maxOrder)
                end
                C = zeros(order+1, order+1)
                S = zeros(order+1, order+1)
            elseif data[1] == "gfc"
                n = parse(Int64, data[2])
                m = parse(Int64, data[3])
                if n ≤ order && m ≤ order
                    C[n+1,m+1] = parse(Float64, data[4])
                    S[n+1,m+1] = parse(Float64, data[5])
                    if m == order && n == order
                        break # skip the rest of the file
                    end
                end
            end
        end
    end
    C[1, 1] = 1.0  # Just in case
    close(file)

    # ================= PRECOMPUTATIONS NEW ALGORITHM ======================= #
    # pre-compute the recursion coefficients corresponding to equations 19 and 22
    # from Holmes and Featherstone paper
    # for cache efficiency, elements are stored in the same order they will be used
    # later on, i.e. from rightmost column to leftmost column
    degree = order
    gnmOj = Float64[]; hnmOj = Float64[]; enm = Float64[]; # For when degree < 2
    for m in degree:-1:0
        j = 2.0 - 1.0*(m > 0)
        for n = max(2, m + 1):degree
            # 	for n = max(1, m + 1):degree
            f = (n - m)*(n + m + 1.0)
            append!(gnmOj, 2.0*(m + 1.0) / sqrt(j*f))
            append!(hnmOj, sqrt((n + m + 2.0)*(n - m - 1.0) / (j*f)))
            append!(enm, sqrt(f / j))
        end
    end
    #gnmOj = gnmOj[2:end]
    #hnmOj = hnmOj[2:end]
    #enm = enm[2:end]

    # scaled sectorial terms corresponding to equation 28 in Holmes and Featherstone paper
    sectorial = zeros(degree + 1)
    sectorial[1] = scalb(1.0, -SCALING)
    sectorial[2] = sqrt(3.0)*sectorial[1]
    for m = 2:degree
        sectorial[m+1] = sqrt((2.0*m + 1.0) / (2.0*m))*sectorial[m]
    end

    return GravityHarmonics(μ, order, R, S, C, SCALING, gnmOj, hnmOj, enm, sectorial,
        zeros(3), ones(degree+1), [1.0; zeros(degree)], zeros(degree+1), zeros(3, 3), zeros(degree+1), zeros(degree+1), zeros(degree+1), zeros(degree+1))
end

@inline scalb(x, n) = x*(2^n)

#  gravity(GH, pos)
#  This function computes the gravitational acceleration of a gravitational
#  potential developed in spherical harmonics. This function is much faster
#  than the equivalent one in STAR, and robust to very high order/degree.
#
#  INPUTS:
#  GH               Gravity Harmonics data structure (to be obtained with GravityModel(...))
#  pos  [m]         Position of the vehicle in planet centered, planet fixed reference frame
#
#  OUTPUTS:
#  acc  [m/s^2]     Gravity acceleration vector in planet fixed reference frame
#
# Reference:
# [1] S. A. Holmes, W. E. Featherstone, A unified approach to the Clenshaw
# summation and the recursive computation of very high degree and order
# normalised associated Legendre functions. Journal of Geodesy (2002) 76:
# 279–299. DOI 10.1007/s00190-002-0216-2.
#
# Adapted from OREKIT
#
function gravity(GH::GravityHarmonics, pos::AbstractVector{T}) where T
    g = zero(pos)
    gravity!(GH, pos, g)
    return g
end

function gravity!(GH::GravityHarmonics, pos::AbstractVector{T}, g::AbstractVector{T}) where T

    # Compute polar coordinates
    x, y, z = pos
    ρ2 = x*x + y*y
    r2 = ρ2 + z*z
    r = sqrt(r2)
    ρ = sqrt(ρ2)
    t = z/r         # cos(theta), where theta is the polar angle
    u = ρ/r         # sin(theta), where theta is the polar angle
    tOu = z/max(ρ, eps(T))

    createDistancePowersArray!(GH, GH.Rref/r)   # compute distance powers
    createCosSinArrays!(GH, x/ρ, y/ρ)           # compute longitude cosines/sines

    # outer summation over order
    index = 1
    U = 0.0
    fill!(GH.dU, 0.0)
    # fill!(GH.pnm0, 0.0)
    # fill!(GH.pnm1, 0.0)
    fill!(GH.pnm0Plus1, 0.0)
    fill!(GH.pnm0Plus2, 0.0)
    pnm0 = GH.pnm0; pnm1 = GH.pnm1; pnm0Plus1 = GH.pnm0Plus1; pnm0Plus2 = GH.pnm0Plus2
    @inbounds for m in GH.order:-1:0

        # compute tesseral terms with derivatives
        index = computeTesseral!(GH, m, index, t, u, tOu, pnm0, pnm1, pnm0Plus1, pnm0Plus2)

        if m ≤ GH.order
            # compute contribution of current order to field (equation 5 of the paper)
            # inner summation over degree, for fixed order
            sumDegreeS = 0.0;        sumDegreeC = 0.0
            dSumDegreeSdR = 0.0;     dSumDegreeCdR = 0.0
            dSumDegreeSdTheta = 0.0; dSumDegreeCdTheta = 0.0
            mp1 = m + 1

            @inbounds for n in max(2, m):GH.order # CAUTION: does this work for 1,1?
                np1 = n + 1
                qSnm = GH.aOrN[np1]*GH.S[np1, mp1]
                qCnm = GH.aOrN[np1]*GH.C[np1, mp1]
                nOr = n/r
                s0 = pnm0[np1]*qSnm
                c0 = pnm0[np1]*qCnm
                s1 = pnm1[np1]*qSnm
                c1 = pnm1[np1]*qCnm
                sumDegreeS += s0
                sumDegreeC += c0
                dSumDegreeSdR -= nOr*s0
                dSumDegreeCdR -= nOr*c0
                dSumDegreeSdTheta += s1
                dSumDegreeCdTheta += c1
            end

            # Contribution to outer summation over order
            sML = GH.sinλ[mp1]
            cML = GH.cosλ[mp1]
            U = muladd(U, u, muladd(sML, sumDegreeS, cML*sumDegreeC))
            GH.dU[1] = GH.dU[1]*u + sML*dSumDegreeSdR + cML*dSumDegreeCdR
            GH.dU[2] = GH.dU[2]*u + sML*dSumDegreeSdTheta + cML*dSumDegreeCdTheta
            GH.dU[3] = GH.dU[3]*u + m*(cML*sumDegreeS - sML*sumDegreeC)
        end

        # Rotate the recursion arrays
        pnm0Plus2, pnm0Plus1, pnm0 = pnm0Plus1, pnm0, pnm0Plus2
    end

    # Scale back
    if GH.SCALING > 0
        U = scalb(U, GH.SCALING)
        GH.dU .= scalb.(GH.dU, GH.SCALING)
    end

    # apply the global mu/r factor
    muOr = GH.μ/r
    GH.dU[1] = GH.dU[1]*muOr - U*muOr/r
    GH.dU[2] *= -muOr
    GH.dU[3] *= muOr

    # Convert Gradient from Spherical to Cartesian Coordinates
    rI = x/r; rJ = y/r
    coeff11 = -GH.μ/(r2*r)*GH.C[1, 1]
    GH.T[1, 1] = rI;   GH.T[1, 2] = -rI*t/ρ;   GH.T[1, 3] = -rJ*r/ρ2
    GH.T[2, 1] = rJ;   GH.T[2, 2] = -rJ*t/ρ;   GH.T[2, 3] = rI*r/ρ2
    GH.T[3, 1] = t;    GH.T[3, 2] = ρ/r2
    mul!(g, GH.T, GH.dU)

    # Add 2-body term (C(1,1) added for compatibility with magnetic field model)
    @inbounds for k in eachindex(g);
        g[k] += coeff11*pos[k];
    end
    return
end

# This function computes normalized associated legendre functions
function computeTesseral!(GH::GravityHarmonics, m, index, t, u, tOu, pnm0, pnm1, pnm0Plus1, pnm0Plus2)
    nmax = max(2, m)

    # initialize recursion from sectorial terms
    n = nmax
    if n == m
        @inbounds pnm0[n+1] = GH.sectorial[n+1]
        n += 1
    end

    # compute tesseral values
    localIndex = index
    u2 = u*u
    @inbounds for k in n+1:GH.order+1
        # value (equation 27 of the paper)
        pnm0[k] = GH.gnmOj[localIndex]*t*pnm0Plus1[k] - GH.hnmOj[localIndex]*u2*pnm0Plus2[k]
        localIndex += 1
    end

    # initialize recursion from sectorial terms
    n = nmax
    if n == m
        @inbounds pnm1[n+1] = m*tOu*pnm0[n+1]
        n += 1
    end

    # compute tesseral values and derivatives with respect to polar angle
    localIndex = index
    @inbounds for k in n+1:GH.order+1
        # first derivative (equation 30 of the paper)
        pnm1[k] = m*tOu*pnm0[k] - GH.enm[localIndex]*u*pnm0Plus1[k]
        localIndex += 1
    end
    return localIndex
end

# This function computes (Rref/r)^n
@inline function createDistancePowersArray!(GH::GravityHarmonics, aOr)
    if GH.order > 0
        GH.aOrN[2] = aOr

        # fill up array
        @inbounds for n in 2:GH.order
            p = fld(n, 2) + 1
            q = n - p + 2
            GH.aOrN[n+1] = GH.aOrN[p]*GH.aOrN[q]
        end
    end
    return
end

# This function computes cos(m*lon) and sin(m*lon)
@inline function createCosSinArrays!(GH::GravityHarmonics, cosLambda, sinLambda)

    # Initialize arrays
    if GH.order > 0
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

            GH.cosλ[m+1] = GH.cosλ[p]*GH.cosλ[q] - GH.sinλ[p]*GH.sinλ[q]
            GH.sinλ[m+1] = GH.sinλ[p]*GH.cosλ[q] + GH.cosλ[p]*GH.sinλ[q]
        end
    end
    return
end
