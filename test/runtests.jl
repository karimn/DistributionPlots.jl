using DistributionPlots
# Tables and Statistics are deps of DistributionPlots but not direct deps of
# test/Project.toml; re-export their bindings into Main here so included test
# files can use them without a separate `using Tables`/`using Statistics`.
using DistributionPlots: Tables, Statistics
using .Statistics: mean, median, quantile
using Test

@testset "DistributionPlots.jl" begin
    @testset "smoke" begin
        @test DistributionPlots isa Module
    end

    include("test_interface.jl")
    include("test_density.jl")
    include("test_intervals.jl")
    include("test_slabs.jl")
end
