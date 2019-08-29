import AstroBase
import JPLEphemeris
using OptionalData: get

export Gravity, UniformGravity, J2Gravity, ThirdBody

abstract type Gravity <: Force end

struct UniformGravity <: Gravity end

function evaluate!(::UniformGravity, δv, t, ep, r, v, params, propagator)
    μ = grav_param(propagator.center)
    rm = norm(r)
    δv .+= -μ .* r ./ rm^3
end

struct J2Gravity <: Gravity end

function evaluate!(::J2Gravity, δv, t, ep, r, v, params, propagator)
    center = propagator.center
    μ = grav_param(center)
    J₂ = j2(center)
    rm = norm(r)
    pj = -3 / 2 * μ * J₂ * mean_radius(center)^2 / (rm^5)

    δv[1:2] .+= -μ .* r[1:2] ./ rm^3 .+ pj .* r[1:2] .*
        (1.0 - 5.0 * r[3]^2 / rm^2)
    δv[3] += -μ * r[3]/rm^3 + pj * r[3] * (3.0 - 5.0 * r[3]^2 / rm^2)
end

struct ThirdBody{T<:AbstractEphemeris} <: Gravity
    ephemeris::T
    bodies::Vector{AstroBase.CelestialBody}
end
ThirdBody(ephemeris::AbstractEphemeris, bodies...) = ThirdBody(ephemeris, collect(bodies))

function evaluate!(tb::ThirdBody, δv, t, ep, r, v, params, propagator)
    rc3 = zeros(3)
    rs3 = zeros(3)
    for body in tb.bodies
        μ = grav_param(body)
        JPLEphemeris.position!(rc3, tb.ephemeris, TDBEpoch(ep), propagator.center, body)
        rs3 .= rc3 .- r
        δv .+= μ * (rs3 ./ norm(rs3)^3 .- rc3 ./ norm(rc3)^3)
        fill!(rc3, 0.0)
    end
end
