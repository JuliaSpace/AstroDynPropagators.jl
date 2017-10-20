using Parameters

import DifferentialEquations: ODEProblem, Vern9, solve, terminate!,
    ContinuousCallback, CallbackSet, set_proposed_dt!
import DifferentialEquations.OrdinaryDiffEq: OrdinaryDiffEqAdaptiveAlgorithm

export ODE

@with_kw struct ODE{F<:Frame,C<:CelestialBody} <: Propagator
    frame::Type{F} = GCRF
    center::Type{C} = Earth
    forces::Vector{Force} = [UniformGravity()]
    minstep::Float64 = 0.0
    maxstep::Float64 = Inf
    algorithm::OrdinaryDiffEqAdaptiveAlgorithm = Vern9()
    events::Vector{Event} = Event[]
end

center(ode::ODE{<:Frame,C}) where {C<:CelestialBody} = C

struct ODEParams
    s0::State
    log::Vector{LogEntry}
end

function propagate(p::ODE{F,C}, s0::State, Δt, points) where {
        F<:Frame, C<:CelestialBody}
    s = State(s0, frame=F, body=C)
    t0 = 0.0
    t1 = seconds(Δt)
    y0 = array(s)
    prm = ODEParams(s, LogEntry[])
    callbacks = ContinuousCallback[]
    for (id, evt) in enumerate(p.events)
        push!(callbacks, ContinuousCallback(
            (t, u, int) -> condition(t, u, int, evt, prm, p),
            (int) -> affect!(int, id, evt, prm, p),
            nothing,
        ))
    end
    prob = ODEProblem(
        (t, y, δy) -> rhs!(t, y, δy, prm, p),
        y0, (t0, t1),
        callback=CallbackSet(callbacks...),
    )
    res = solve(prob, p.algorithm;
        dtmin=p.minstep, dtmax=p.maxstep,
        save_everystep=points != :none,
    )
    if res.retcode == :MaxIters
        error("Maximum number of iterations reached. Propagation aborted.")
    elseif res.retcode != :Success
        error("Solver returned error: $(res.retcode)")
    end
    ep1 = epoch(s0) + res.t[end] * seconds
    s1 = State(ep1, res.u[end][1:3], res.u[end][4:6], F, C)
    Trajectory(s, s1, res.t, res.u, prm.log)
end

function rhs!(t, y, δy, params, propagator)
    ep = epoch(params.s0) + t * seconds
    δy[1:3] .= y[4:6]
    @views begin
        r = y[1:3]
        v = y[4:6]
        δv = δy[4:6]
    end
    fill!(δv, 0.0)
    for force in propagator.forces
        evaluate!(force, δv, t, ep, r, v, params, propagator)
    end
    if isrotating(propagator.frame)
        rotational!(δv, t, ep, r, v, params, propagator)
    end
end

function condition(t, y, integrator, evt, params, propagator)
    detect(evt.detector, t, y, params, propagator)
end

function affect!(integrator, idx, evt, params, propagator)
    if !evt.detect_all
        undetected = count_id(idx, params.log) == 0
        !undetected && return
    end

    ep = epoch(params.s0) + integrator.t * seconds
    name = Base.datatype_name(typeof(evt.detector))
    push!(params.log, LogEntry(idx, name, integrator.t, ep))

    if !isnull(evt.updater)
        update!(get(evt.updater), integrator, idx, params, propagator)
        set_proposed_dt!(integrator, 1.0)
    end
end
