export UniformGravity, J2Gravity, ThirdBody

struct UniformGravity{C<:CelestialBody} <: Force end

UniformGravity(::Type{C}) where {C<:CelestialBody} = UniformGravity{C}()

function evaluate!(::UniformGravity{C}, δv, t, r, v) where C<:CelestialBody
    μ = mu(C)
    rm = norm(r)
    δv .+= -μ .* r ./ rm^3
end

struct J2Gravity{C<:CelestialBody} <: Force end

J2Gravity(::Type{C}) where {C<:CelestialBody} = J2Gravity{C}()

function evaluate!(::J2Gravity{C}, δv, t, r, v) where C<:CelestialBody
    μ = mu(C)
    J₂ = j2(C)
    rm = norm(r)
    pj = -3 / 2 * μ * J₂ * mean_radius(C)^2 / (rm^5)

    δv[1:2] .+= -μ .* r[1:2] ./ rm^3 .+ pj .* r[1:2] .* (1.0 - 5.0 * r[3]^2 / rm^2)
    δv[3] += -μ * r[3]/rm^3 + pj * r[3] * (3.0 - 5.0 * r[3]^2 / rm^2)
end

struct ThirdBody <: Force
    bodies::Vector{DataType}
end
