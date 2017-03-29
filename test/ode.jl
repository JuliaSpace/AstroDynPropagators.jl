using AstroDynBase
using AstroDynCoordinates
using JPLEphemeris
using RemoteFiles

de430 = @RemoteFile "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp"
download(de430)
load_ephemeris!(SPK, path(de430))

@testset "ODE Propagator" begin
    o = ODE()
    r0 = [1131.340, -2282.343, 6672.423]km
    v0 = [-5.64305, 4.30333, 2.42879]kps
    Δt = 40*60*seconds
    ep = TTEpoch(2000, 1, 1)
    s0 = State(ep, r0, v0)
    r1 = [-4219.7527, 4363.0292, -3958.7666]km
    v1 = [3.689866, -1.916735, -6.112511]kps
    s1 = State(ep + Δt, r1, v1)
    ode = ODE(maxstep=10)
    tra = propagate(ode, s0, Δt)
    @test isapprox(radius(final(tra)), radius(s1), rtol=1e-4)
    @test isapprox(velocity(final(tra)), velocity(s1), rtol=1e-4)
    
    # s2 = state(s0, Δt/4, ode)
    # s1 = state(s0, Δt, ode)
    # @test s1.rv ≈ [r1; v1]
    # s1 = state(s0, Δe, ode)
    # @test s1.rv ≈ [r1; v1]
    # tra = trajectory(s0, Δe, ode)
    # s1 = tra.trajectory[ep+Δe]
    # @test s1.rv ≈ [r1; v1]
    # @test tra.trajectory[Δt/4] ≈ s2
    # ode = ODE(discontinuities=[Discontinuity(PericenterEvent(), Stop())], maxstep=100)
    # s1 = state(s0, period(s0), ode)
    # ano = trueanomaly(s1)
    # ano = ano > π ? abs(ano - 2π) : ano
    # @test round(ano, 10) ≈ 0.0
    # ode = ODE(discontinuities=[Discontinuity(ApocenterEvent(), Stop())], maxstep=100)
    # s1 = state(s0, period(s0), ode)
    # ano = trueanomaly(s1)
    # @test ano ≈ π
    # ode = ODE(discontinuities=[Discontinuity(ApocenterEvent(), Abort())], maxstep=100)
    # @test_throws PropagatorAbort state(s0, period(s0), ode)
    # ode = ODE(discontinuities=[Discontinuity(StartEvent(), ImpulsiveManeuver(along=-1))])
    # @test_throws PropagatorAbort state(s0, period(s0), ode)
    # ode = ODE(events=[ApocenterEvent()], maxstep=10)
    # tra = trajectory(s0, period(s0)*4, ode)
    # for t in tra.events[:index]
    #     ano = trueanomaly(tra.trajectory[t])
    #     @test ano ≈ π
    # end
end
