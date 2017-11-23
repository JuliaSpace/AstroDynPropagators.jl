using AstroDynBase
using Parameters

import AstroDynBase: epoch
import OrdinaryDiffEq: terminate!

export Detector, Updater, Event, detect, update!,
    Apocenter, Pericenter, Timed, Impact, Height,
    Abort, Stop, ImpulsiveManeuver

abstract type Detector end
abstract type Updater end

@with_kw struct Event
    detector::Detector
    updater::Nullable{Updater} = Nullable{Updater}()
    detect_all::Bool = isnull(updater)
end

function detect(det::Detector, t, y, params, propagator)
    error("No 'detect' method defined for '$(Base.datatype_name(typeof(det)))'")
end

struct Timed{T<:Number} <: Detector
    time::T
end

function detect(timed::Timed, t, y, params, propagator)
    t - timed.time
end

struct Apocenter <: Detector end

function detect(::Apocenter, t, y, params, propagator)
    el = keplerian(y, μ(center(propagator)))
    ano = el[6]
    if ano > pi/2
        ano = abs(ano - pi)
    elseif ano < -pi/2
        ano = -abs(ano + pi)
    end
    isretrograde(el[3]) ? ano : -ano
end

struct Pericenter <: Detector end

function detect(::Pericenter, t, y, params, propagator)
    el = keplerian(y, μ(center(propagator)))
    isretrograde(el[3]) ? -el[6] : el[6]
end

struct Impact <: Detector end

function detect(::Impact, t, y, params, propagator)
    -(norm(y[1:3]) - mean_radius(center(propagator)))
end

struct Height{T<:Number} <: Detector
    height::T
    ascending::Bool
end

function detect(h::Height, t, y, params, propagator)
    height = norm(y[1:3]) - mean_radius(center(propagator)) - h.height
    h.ascending ? height : -height
end

@with_kw struct Abort <: Updater
    msg::String = "Propagation aborted."
    num::Int = 1
end

function update!(u::Abort, integrator, id, params, propagator)
    if count_id(id, params.log) == u.num
        error(u.msg)
    end
end

@with_kw struct Stop <: Updater
    num::Int = 1
end

function update!(u::Stop, integrator, id, params, propagator)
    if count_id(id, params.log) == u.num
        terminate!(integrator)
    end
end

struct ImpulsiveManeuver{T<:Number} <: Updater
    Δv::Array{T}
end

function ImpulsiveManeuver(;radial=0.0, along=0.0, cross=0.0)
    ImpulsiveManeuver([radial, along, cross])
end

function update!(man::ImpulsiveManeuver, integrator, id, params, propagator)
    rot = Rotation(RAC, propagator.frame, integrator.u[1:3], integrator.u[4:6])
    Δv, _ = rot(man.Δv, zeros(3))
    integrator.u[4:6] .+= Δv
end
