export UniformGravity, J2Gravity

struct UniformGravity{C<:CelestialBody} <: Force end

UniformGravity(::Type{C}) where {C<:CelestialBody} = UniformGravity{C}()

function evaluate!(::UniformGravity{C}, f, y) where C<:CelestialBody
    r = norm(y[1:3])
    f[1:3] .+= y[4:6]
    f[4:6] .+= -μ(C) .* y[1:3] ./ r^3
end

struct J2Gravity{C<:CelestialBody} <: Force end

J2Gravity(::Type{C}) where {C<:CelestialBody} = J2Gravity{C}()

function evaluate!(::J2Gravity{C}, f, y) where C<:CelestialBody
    r = norm(y[1:3])
    pj = -3/2 * μ(C) * j2(C) * mean_radius(C)^2 / (r^5)

    f[1:3] .+= y[4:6]
    f[4:5] .+= -μ(C) .* y[1:2] ./ r^3 .+ pj .* y[1:2] .* (1.0 - 5.0 * z^2 / r^2)
    f[6] += -μ(C) * y[3]/r^3 + pj * y[3] * (3.0 - 5.0 * z^2 / r^2)
end
