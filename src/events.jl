using AstroDynBase
using Parameters

import AstroDynBase: epoch

abstract type Detector end
abstract type Updater end

@with_kw struct Event
    detector::Detector
    updater::Nullable{Updater} = nothing
    multi::Bool = isnull(updater)
end

struct Timer{T<:Number} <: Detector
    time::T
end
struct Start <: Detector end
struct End <: Detector end
struct Apocenter <: Detector end
struct Pericenter <: Detector end
struct Impact <: Detector end

struct Abort <: Updater
    num::Int
    msg::String
    Abort(num=1, msg="Propagation aborted.") = Abort(num, msg)
end

struct LogEntry
    id::Int
    detector::Symbol
    t::Epoch
end

id(l::LogEntry) = l.id
detector(l::LogEntry) = l.detector
epoch(l::LogEntry) = l.epoch
