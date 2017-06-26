using AstroDynBase
using Parameters

import AstroDynBase: epoch
import DifferentialEquations: terminate!

export Detector, Updater, Event, detect, update!,
    Apocenter, Pericenter, Timed,
    Abort, Stop

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
    -norm(y[1:3]) - mean_radius(center(propagator))
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

struct LogEntry
    id::Int
    detector::Symbol
    time::Float64
    ep::Epoch
end

id(l::LogEntry) = l.id
detector(l::LogEntry) = l.detector
epoch(l::LogEntry) = l.epoch
count_id(idx, log) = count(x->id(x) == idx, log)
