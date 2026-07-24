using DistributionPlots
using Test

@testset "DistributionPlots.jl" begin
    @testset "smoke" begin
        @test DistributionPlots isa Module
    end
end
