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
        solout=solout!,
        points=:all,
        maxstep=p.maxstep,
        numstep=p.numstep,
    )
    ep1 = epoch(s0) + tout[end] * seconds
    s1 = State(ep1, yout[end][1:3]km, yout[end][4:6]kps, F, C)
    x = map(v -> v[1]km, yout)
    y = map(v -> v[2]km, yout)
    z = map(v -> v[3]km, yout)
    vx = map(v -> v[4]kps, yout)
    vy = map(v -> v[5]kps, yout)
    vz = map(v -> v[6]kps, yout)
    Trajectory(s0, s1, tout * seconds, x, y, z, vx, vy, vz)
end

function rhs!(f, t, y, params, propagator)
    fill!(f, 0.0)
    for force in propagator.forces
        evaluate!(force, f, t, y)
    end
end

function solout!(told, t, y, contd, params, propagator)
    firststep = told ≈ t
    if (i, evt) in enumerate(params.events)
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
