using AstroDynBase
using AstroDynCoordinates

@testset "Kepler Propagator" begin
    r0 = [1131.340, -2282.343, 6672.423]
    v0 = [-5.64305, 4.30333, 2.42879]
    Δt = 40*60seconds
    ep = TTEpoch(2017, 3, 29)
    ep1 = ep + Δt
    s0 = State(ep, r0, v0)
    r1 = [-4219.7527, 4363.0292, -3958.7666]
    v1 = [3.689866, -1.916735, -6.112511]
    s1 = State(ep1, r1, v1)
    k = Kepler()
    tra = propagate(k, s0, Δt, points=:none)
    @test isapprox(radius(final(tra)), radius(s1), rtol=1e-4)
    @test isapprox(velocity(final(tra)), velocity(s1), rtol=1e-4)

    s2 = final(propagate(k, s0, Δt/4))
    tra = propagate(k, s0, Δt)
    @test final(tra) == tra(Δt)
    @test final(tra) == tra(ep1)
    s3 = tra(Δt)
    @test isapprox(radius(s3), radius(s1), rtol=1e-4)
    @test isapprox(velocity(s3), velocity(s1), rtol=1e-4)
    s4 = tra(Δt/4)
    @test isapprox(radius(s4), radius(s2), rtol=1e-4)
    @test isapprox(velocity(s4), velocity(s2), rtol=1e-4)
    @test tra[1] == initial(tra)
    @test tra[k.points] == final(tra)
    @test tra[end] == final(tra)
    tra = propagate(Kepler(), s0)
    @test tra[end] ≈ State(ep + period(s0), r0, v0)
end
