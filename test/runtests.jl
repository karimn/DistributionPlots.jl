using DistributionPlots
using Test

@testset "DistributionPlots.jl" begin
    @testset "smoke" begin
        @test DistributionPlots isa Module
    end

    include("test_interface.jl")
    include("test_density.jl")
end
