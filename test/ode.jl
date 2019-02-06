using AstroBase: sun, moon
using AstroDynBase
using AstroDynCoordinates
using JPLEphemeris
using RemoteFiles

de430 = @RemoteFile("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/" *
    "planets/de430.bsp")
download(de430)
load_ephemeris!(SPK, path(de430))

@testset "ODE Propagator" begin
    r0 = [1131.340, -2282.343, 6672.423]
    v0 = [-5.64305, 4.30333, 2.42879]
    ep = TTEpoch(2000, 1, 1)
    s0 = State(ep, r0, v0)
    s1 = State(ep + period(s0), r0, v0)
    ode = ODE(maxstep=100.0)
    tra = propagate(ode, s0)
    @test final(tra) ≈ s1

    s0_rot = State(s0, frame=IAUEarth)
    s1_rot = State(State(ep + period(s0), r0, v0), frame=IAUEarth)
    ode = ODE(frame=IAUEarth, maxstep=100.0)
    tra = propagate(ode, s0)
    @test initial(tra) == s0_rot
    @test final(tra) ≈ s1_rot

    # Reference values from Orekit
    r1 = [-4255.223590627231, 4384.471704756651, -3936.1350079623207]
    v1 = [3.6559899898490054, -1.884445831960271, -6.123308149589636]
    Δt = 40*60*seconds
    s1 = State(ep + Δt, r1, v1)
    ode = ODE(
        forces=[J2Gravity()],
        maxstep=100.0,
    )
    tra = propagate(ode, s0, Δt)
    @test final(tra) ≈ s1

    # Reference values from Orekit
    r1 = [-4219.7545636149, 4363.0305735489, -3958.767123328]
    v1 = [3.689863232, -1.9167326384, -6.1125111597]
    ep = TDBEpoch(2020, 1, 1)
    s0 = State(ep, r0, v0)
    s1 = State(ep + Δt, r1, v1)
    ode = ODE(
        forces=[UniformGravity(), ThirdBody(sun, moon)],
        maxstep=100.0,
    )
    tra = propagate(ode, s0, Δt)
    @test final(tra) ≈ s1
end
