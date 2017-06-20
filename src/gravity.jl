export Gravity, UniformGravity, J2Gravity, ThirdBody

abstract type Gravity <: Force end

struct UniformGravity <: Gravity end

function evaluate!(::UniformGravity, δv, t, ep, r, v, params, propagator)
    μ = mu(propagator.center)
    rm = norm(r)
    δv .+= -μ .* r ./ rm^3
end

struct J2Gravity <: Gravity end

function evaluate!(::J2Gravity, δv, t, ep, r, v, params, propagator)
    center = propagator.center
    μ = mu(center)
    J₂ = j2(center)
    rm = norm(r)
    pj = -3 / 2 * μ * J₂ * mean_radius(center)^2 / (rm^5)

    δv[1:2] .+= -μ .* r[1:2] ./ rm^3 .+ pj .* r[1:2] .*
        (1.0 - 5.0 * r[3]^2 / rm^2)
    δv[3] += -μ * r[3]/rm^3 + pj * r[3] * (3.0 - 5.0 * r[3]^2 / rm^2)
end

struct ThirdBody <: Gravity
    bodies::Vector{DataType}
end
ThirdBody(bodies...) = ThirdBody(collect(bodies))

function evaluate!(tb::ThirdBody, δv, t, ep, r, v, params, propagator)
    rc3 = zeros(3)
    rs3 = zeros(3)
    for body in tb.bodies
        μ = mu(body)
        rc3 .= position(ep, propagator.center, body)
        rs3 .= rc3 .- r
        δv .+= μ * (rs3 ./ norm(rs3)^3 .- rc3 ./ norm(rc3)^3)
    end
end
