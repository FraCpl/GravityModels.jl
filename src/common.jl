function gravity(GH::G, pos::AbstractVector{T}) where {T, G<:AbstractGravity}
    g = zero(pos)
    x, y, z = pos
    g[1], g[2], g[3] = gravity(GH, x, y, z)
    return g
end

function gravity!(GH::G, pos::AbstractVector{T}, g::AbstractVector{T}) where {T, G<:AbstractGravity}
    x, y, z = pos
    g[1], g[2], g[3] = gravity(GH, x, y, z)
    return g
end

function gravity!(GH::GravityHarmonics, R_PX::AbstractMatrix{M}, pos_X::AbstractVector{T}, g_X::AbstractVector{T}) where {T, M}
    # Compute position coordinates in planet-centric frame P
    X, Y, Z = pos_X
    x = R_PX[1, 1] * X + R_PX[1, 2] * Y + R_PX[1, 3] * Z
    y = R_PX[2, 1] * X + R_PX[2, 2] * Y + R_PX[2, 3] * Z
    z = R_PX[3, 1] * X + R_PX[3, 2] * Y + R_PX[3, 3] * Z

    # Compute gravity
    gx, gy, gz = gravity(GH, x, y, z)

    # Re-project gravity into input frame X
    g_X[1] = R_PX[1, 1] * gx + R_PX[2, 1] * gy + R_PX[3, 1] * gz
    g_X[2] = R_PX[1, 2] * gx + R_PX[2, 2] * gy + R_PX[3, 2] * gz
    g_X[3] = R_PX[1, 3] * gx + R_PX[2, 3] * gy + R_PX[3, 3] * gz
    return g_X
end

gravity!(GH::GravityKeplerian, R_PX::AbstractMatrix{M}, pos_X::AbstractVector{T}, g_X::AbstractVector{T}) where {T, M} = gravity!(GH, pos_X, g_X)

function gravityThirdBody!(GH::GravityKeplerian, R_PX::AbstractMatrix{M}, pos_X::AbstractVector{T}, posBody_X::AbstractVector{T}, g_X::AbstractVector{T}) where {T, M}
    x, y, z = pos_X
    X, Y, Z = posBody_X
    dx = x - X
    dy = y - Y
    dz = z - Z

    # Gravity of the third body on the spacecraft
    gx, gy, gz = gravity(GH, dx, dy, dz)

    # Gravity of the third body on the frame origin X
    gX, gY, gZ = gravity(GH, -X, -Y, -Z)

    g_X[1] = gx - gX
    g_X[2] = gy - gY
    g_X[3] = gz - gZ
    return g_X
end

function gravityThirdBody!(GH::GravityHarmonics, R_PX::AbstractMatrix{M}, pos_X::AbstractVector{T}, posBody_X::AbstractVector{T}, g_X::AbstractVector{T}) where {T, M}
    x0, y0, z0 = pos_X
    X0, Y0, Z0 = posBody_X

    # Project coordinates in P
    x = R_PX[1, 1] * x0 + R_PX[1, 2] * y0 + R_PX[1, 3] * z0
    y = R_PX[2, 1] * x0 + R_PX[2, 2] * y0 + R_PX[2, 3] * z0
    z = R_PX[3, 1] * x0 + R_PX[3, 2] * y0 + R_PX[3, 3] * z0

    X = R_PX[1, 1] * X0 + R_PX[1, 2] * Y0 + R_PX[1, 3] * Z0
    Y = R_PX[2, 1] * X0 + R_PX[2, 2] * Y0 + R_PX[2, 3] * Z0
    Z = R_PX[3, 1] * X0 + R_PX[3, 2] * Y0 + R_PX[3, 3] * Z0

    dx = x - X
    dy = y - Y
    dz = z - Z

    # Gravity of the third body on the spacecraft
    gx, gy, gz = gravity(GH, dx, dy, dz)

    # Gravity of the third body on the frame origin X
    gX, gY, gZ = gravity(GH, -X, -Y, -Z)

    dgx = gx - gX
    dgy = gy - gY
    dgz = gz - gZ

    # Reproject back to X
    g_X[1] = R_PX[1, 1] * dgx + R_PX[2, 1] * dgy + R_PX[3, 1] * dgz
    g_X[2] = R_PX[1, 2] * dgx + R_PX[2, 2] * dgy + R_PX[3, 2] * dgz
    g_X[3] = R_PX[1, 3] * dgx + R_PX[2, 3] * dgy + R_PX[3, 3] * dgz
    return g_X
end
