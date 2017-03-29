export UniformGravity, J2Gravity

struct UniformGravity{C<:CelestialBody} <: Force end

UniformGravity(::Type{C}) where {C<:CelestialBody} = UniformGravity{C}()

function evaluate!(::UniformGravity{C}, f, t, y) where C<:CelestialBody
    μ = ustrip(mu(C))
    r = norm(y[1:3])
    f[1:3] .+= y[4:6]
    f[4:6] .+= -μ .* y[1:3] ./ r^3
end

struct J2Gravity{C<:CelestialBody} <: Force end

J2Gravity(::Type{C}) where {C<:CelestialBody} = J2Gravity{C}()

function evaluate!(::J2Gravity{C}, f, t, y) where C<:CelestialBody
    μ = ustrip(mu(C))
    J₂ = j2(C)
    rm = ustrip(mean_radius(C))
    r = norm(y[1:3])
    pj = -3/2 * μ * J₂ * mean_radius(C)^2 / (r^5)

    f[1:3] .+= y[4:6]
    f[4:5] .+= -μ .* y[1:2] ./ r^3 .+ pj .* y[1:2] .* (1.0 - 5.0 * z^2 / r^2)
    f[6] += -μ * y[3]/r^3 + pj * y[3] * (3.0 - 5.0 * z^2 / r^2)
end
