using AstroDynBase
using AstroDynCoordinates

export Propagator, propagate

import Base: show

abstract type Propagator end
abstract type Force end

function propagate(p::Propagator, s0::AbstractState, Δt; points=:all)
    Δt = in_seconds.(Δt)
    propagate(p, s0, Δt, points)
end

function propagate(p::Propagator, s0::State; points=:all)
    Δt = period(s0)
    propagate(p, s0, Δt, points)
end
