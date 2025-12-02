# Geodetic coordinates to ECEF geocentric coordinates
#
# Reference:
# [1] Zhu, Conversion of Earth-Cented Earth-Fixed Coordinates o Geodetic
#     Coordinates, IEEE Transactions on Aerospace and Electronic Systems,
#     Vol. 30, No. 3, July 1994.
#     https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=303772
#
# Author: F. Capolupo
# European Space Agency, 2021
function geodetic2geocentric(lat, lon, h, a = 6378137.0, f = 298.257223563)
    pos = Vector{typeof(lat)}(undef, 3)
    geodetic2geocentric!(pos, lat, lon, h, a, f)
    return pos
end

function geodetic2geocentric!(pos, lat, lon, h, a = 6378137.0, f = 298.257223563)
    fi = 1/f
    e2 = fi*(2 - fi)

    sφ, cφ = sincos(lat)
    sλ, cλ = sincos(lon)

    Rc = a/sqrt(1 - e2*sφ*sφ)
    p = (Rc + h)*cφ
    pos[1] = p*cλ
    pos[2] = p*sλ
    pos[3] = ((1 - e2)*Rc + h)*sφ
    return
end


function gravityMap(GH::GravityHarmonics, altitude = 100e3; N = 400)
    r = GH.Rref + altitude
    lat = range(-90, 90, N)[2:(end-1)]
    lon = range(-180, 180, 2N)[1:(end-1)]
    return lat,
    lon,
    [
        gravity(GH, r*[cosd(latk)*cosd(lonk); cosd(latk)*sind(lonk); sind(latk)]) for
        lonk in lon, latk in lat
    ]
end
