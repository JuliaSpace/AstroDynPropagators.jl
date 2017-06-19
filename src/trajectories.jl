using AstroDynBase
using AstroDynCoordinates
using SmoothingSplines

import Base: getindex, endof, show

export Trajectory, initial, final, state, events, times

abstract type AbstractTrajectory end

struct Trajectory <: AbstractTrajectory
    initial::AbstractState
    final::AbstractState
    times::Vector
    splines::Vector{SmoothingSpline}
    events::Vector{LogEntry}
end

initial(tra::AbstractTrajectory) = tra.initial
final(tra::AbstractTrajectory) = tra.final
events(tra::AbstractTrajectory) = tra.events
times(tra::AbstractTrajectory) = tra.times

function Trajectory(initial, final, events=LogEntry[])
    Trajectory(initial, final, SmoothingSpline[], events)
end

function Trajectory(initial, final, t, x, y, z, vx, vy, vz, events=LogEntry[])
    splines = Array{SmoothingSpline}(6)
    splines[1] = fit(SmoothingSpline, t, x, 0.0)
    splines[2] = fit(SmoothingSpline, t, y, 0.0)
    splines[3] = fit(SmoothingSpline, t, z, 0.0)
    splines[4] = fit(SmoothingSpline, t, vx, 0.0)
    splines[5] = fit(SmoothingSpline, t, vy, 0.0)
    splines[6] = fit(SmoothingSpline, t, vz, 0.0)
    Trajectory(initial, final, t, splines, events)
end

function show(io::IO, tra::Trajectory)
    println(io, "Trajectory")
    println(io, " Start date: $(epoch(initial(tra)))")
    println(io, " End date:   $(epoch(final(tra)))")
end

interpolate(tra::Trajectory, time) = SmoothingSplines.predict.(tra.splines, time)

function state(tra::Trajectory, time)
    rv = interpolate(tra, time)
    r = rv[1:3] * unit(tra.x[1])
    v = rv[4:6] * unit(tra.vx[1])
    ep1 = epoch(initial(tra)) + second(time)
    f = frame(initial(tra))
    b = body(initial(tra))
    return State(ep1, r, v, f, b)
end

function state(tra::Trajectory, ep::Epoch)
    ep0 = epoch(initial(tra))
    time = typeof(ep0)(ep) - ep0
    tra(second(time))
end

(tra::Trajectory)(time) = state(tra, time)

getindex(tra::Trajectory, idx) = state(tra, tra.times[idx])
endof(tra::Trajectory) = length(tra.times)
