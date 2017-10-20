@testset "Events" begin
    r = [
    6068279.27,
    -1692843.94,
    -2516619.18,
    ]/1000

    v = [
    -660.415582,
    5495.938726,
    -5303.093233,
    ]/1000

    t = UTCEpoch("2016-05-30T12:00:00.000")
    iss = State(t, r, v)

    ode = ODE(events=[Event(detector=Pericenter())])
    tra = propagate(ode, iss)
    @test AstroDynPropagators.count_id(1, tra.events) == 1

    ode = ODE(events=[Event(detector=Pericenter())])
    tra = propagate(ode, iss, period(iss) * 2)
    @test AstroDynPropagators.count_id(1, tra.events) == 2

    ode = ODE(events=[Event(detector=Pericenter(), updater=Stop())])
    tra = propagate(ode, iss)
    @test isapprox(trueanomaly(final(tra)), 0.0, atol=1e-12)

    ode = ODE(events=[Event(detector=Apocenter())])
    tra = propagate(ode, iss)
    @test AstroDynPropagators.count_id(1, tra.events) == 1

    ode = ODE(events=[Event(detector=Apocenter())])
    tra = propagate(ode, iss, period(iss) * 2)
    @test AstroDynPropagators.count_id(1, tra.events) == 2

    ode = ODE(events=[Event(detector=Apocenter(), updater=Stop())])
    tra = propagate(ode, iss)
    @test abs(trueanomaly(final(tra))) ≈ π

    ode = ODE(events=[Event(detector=Timed(60.0))])
    tra = propagate(ode, iss)
    @test AstroDynPropagators.count_id(1, tra.events) == 1

    ode = ODE(events=[Event(detector=Timed(60.0), updater=Stop())])
    tra = propagate(ode, iss)
    res = seconds(epoch(final(tra)) - epoch(initial(tra)))
    @test res ≈ seconds(60seconds)

    iss1 = State(t, r, v * 0.7)
    ode = ODE(events=[Event(detector=Height(200.0, false), updater=Stop())])
    tra = propagate(ode, iss1)
    @test norm(position(final(tra))) - mean_radius(body(iss1)) ≈ 200.0

    ode = ODE(events=[Event(detector=Impact(), updater=Stop())])
    tra = propagate(ode, iss1)
    height = norm(position(final(tra))) - mean_radius(body(iss1))
    @test isapprox(height, 0.0, atol=1e-10)

    events = [
        Event(detector=Pericenter(), updater=ImpulsiveManeuver(along=-1.0)),
        Event(detector=Impact(), updater=Abort()),
    ]
    ode = ODE(events=events)
    @test_throws ErrorException propagate(ode, iss)
end
