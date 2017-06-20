using AstroDynPropagators
using Base.Test

@testset "AstroDynPropagators" begin
    include("kepler.jl")
    include("ode.jl")
end
