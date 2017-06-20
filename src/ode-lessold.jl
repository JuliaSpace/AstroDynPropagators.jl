using AstroDynBase
using Dopri

export ODE

struct ODE{F<:Frame, C<:CelestialBody} <: Propagator
    forces::Array{Force}
    events::Array{Event}
    maxstep::Float64
    numstep::Int
end

function ODE(forces...;
    frame::Type{F}=GCRF,
    center::Type{C}=Earth,
    maxstep=0.0,
    numstep=100_000,
    events=Event[],
) where {F<:Frame, C<:CelestialBody}
    if isempty(forces)
        forces = [UniformGravity(C)]
    end
    ODE{F, C}(forces, events, maxstep, numstep)
end

mutable struct ODEParams
    events::Array{Event}
end

function propagate(p::ODE{F,C}, s0::State, Δt, points) where {F<:Frame, C<:CelestialBody}
    s = State(s0, frame=F, body=C)
    t0 = 0.0
    t1 = ustrip(Δt)
    y = array(s)
    tout, yout = dop853((f, t, y) -> rhs!(f, t, y, 0.0, p), y, [t0, t1],
        # solout=solout!,
        points=:all,
        maxstep=p.maxstep,
        numstep=p.numstep,
    )
    ep1 = epoch(s0) + tout[end]
    s1 = State(ep1, yout[end][1:3], yout[end][4:6], F, C)
    x = map(v -> v[1], yout)
    y = map(v -> v[2], yout)
    z = map(v -> v[3], yout)
    vx = map(v -> v[4], yout)
    vy = map(v -> v[5], yout)
    vz = map(v -> v[6], yout)
    Trajectory(s0, s1, tout, x, y, z, vx, vy, vz)
end

function rhs!(f, t, y, params, propagator)
    δv = fill(0.0, 3)
    for force in propagator.forces
        evaluate!(force, δv, t, y[1:3], y[4:6])
    end
    f[1:3] = y[4:6]
    f[4:6] = ustrip(δv)
end

function solout!(told, t, y, contd, params, propagator)
    firststep = told ≈ t
    for (i, evt) in enumerate(params.events)
        isdone(evt) && continue

        if isdetected(evt) && !firststep
        end
    end
end

function state(contd, t, rng)
    y = Array{Float64}(length(rng))
    for i in rng
        y[i] = contd(i, t)
    end
    return y
end
