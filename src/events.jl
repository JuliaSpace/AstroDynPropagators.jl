abstract type Detector end
abstract type Updater end

struct Event
    time::Nullable{Float64}
    detector::Nullable{Detector}
    updater::Nullable{Updater}
    done::Bool
    log::Vector{Float64}
end
Event(detector::Detector, updater=nothing) = Event(nothing, detector, updater, false, [])
Event(time::Float64, updater=nothing) = Event(time, nothing, updater, false, [])
isdone(evt::Event) = evt.done
istimed(evt::Event) = !isnull(evt.time)
isdetected(evt::Event) = !isnull(evt.detector)

time(evt::Event) = !isnull(evt.time) ? get(evt.time) : time(get(evt.detector))

struct Start <: Detector end
struct End <: Detector end
struct Apocenter <: Detector end
struct Pericenter <: Detector end
struct Impact <: Detector end

struct Abort <: Updater
    num::Int
    msg::String
    Abort(num=1, msg="Propagation abortedd") = Abort(num, msg)
end
