using AstroDynModels
using Base.Test

@testset "AstroDynModels" begin
    include("kepler.jl")
    include("ode.jl")
end
