using AstroDynBase
using AstroDynCoordinates
using Parameters
using Unitful

export Kepler

@with_kw struct Kepler <: Propagator
    iterations::Int = 50
    points::Int = 100
    rtol::Float64 = sqrt(eps())
end

show(io::IO, ::Type{Kepler}) = print(io, "Kepler")

function propagate(p::Kepler, s0::State, Δt, points)
    ep1 = epoch(s0) + Δt
    if points == :none
        r1, v1 = kepler(μ(body(s0)), radius(s0), velocity(s0), Δt,
            p.iterations, p.rtol)
        s1 = State(ep1, r1, v1, frame(s0), body(s0))
        Trajectory(s0, s1)
    else
        times = linspace(zero(Δt), Δt, p.points)
        runit = typeof(1.0 * unit(radius(s0)[1]))
        vunit = typeof(1.0 * unit(velocity(s0)[1]))
        x = zeros(runit, p.points)
        y = zeros(runit, p.points)
        z = zeros(runit, p.points)
        vx = zeros(vunit, p.points)
        vy = zeros(vunit, p.points)
        vz = zeros(vunit, p.points)
        for (i, t) in enumerate(times)
            r, v = kepler(μ(body(s0)), radius(s0), velocity(s0), t,
                p.iterations, p.rtol)
            x[i] = r[1]
            y[i] = r[2]
            z[i] = r[3]
            vx[i] = v[1]
            vy[i] = v[2]
            vz[i] = v[3]
        end
        s1 = State(ep1, [x[end], y[end], z[end]], [vx[end], vy[end], vz[end]],
            frame(s0), body(s0))
        Trajectory(s0, s1, times, x, y, z, vx, vy, vz)
    end
end