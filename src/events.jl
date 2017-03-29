abstract type Detector end
abstract type Updater end

struct Event
    time::Nullable{Epoch}
    detector::Nullable{Detector}
    updater::Nullable{Updater}
end
Event(detector::Detector, updater=nothing) = Event(nothing, detector, updater)
Event(time::Epoch, updater=nothing) = Event(time, nothing, updater)

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
