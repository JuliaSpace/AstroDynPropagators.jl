using AstroDynPropagators
using Test

@testset "AstroDynPropagators" begin
    include("kepler.jl")
    include("ode.jl")
    include("events.jl")
end
