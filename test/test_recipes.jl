using DistributionPlots
using DistributionPlots: SlabInterval, _to_dist_args, AnalyticDist, SampleDist
using DistributionPlots: Makie
using .Makie
using Distributions
using Test

@testset "recipe scaffolding" begin
    # shared converter: a single distribution → one position, one dist
    pos, dists = _to_dist_args(Normal(0, 1))
    @test pos == [1.0]
    @test dists[1] isa AnalyticDist

    # a plain sample vector → one SampleDist at position 1
    pos2, dists2 = _to_dist_args(randn(500))
    @test pos2 == [1.0]
    @test dists2[1] isa SampleDist

    # a vector of distributions → positions 1:k
    pos3, dists3 = _to_dist_args([Normal(0, 1), Normal(1, 2), Normal(2, 3)])
    @test pos3 == [1.0, 2.0, 3.0]
    @test all(d isa AnalyticDist for d in dists3)

    # the recipe type exists and carries our attributes with ggdist-style names
    @test SlabInterval <: Makie.AbstractPlot
    attrs = Makie.default_theme(nothing, SlabInterval)
    for a in (:slab_type, :interval, :point, :side, :justification, :scale,
              :normalize, :show_slab, :show_interval, :show_point, :trim, :colormap)
        @test haskey(attrs, a)
    end

    # convert_arguments dispatches through the shared converter
    ca1 = Makie.convert_arguments(SlabInterval, Normal(0, 1))
    @test ca1 == (_to_dist_args(Normal(0, 1)),)

    ca2 = Makie.convert_arguments(SlabInterval, randn(100))
    posr, distsr = ca2[1]
    @test posr == [1.0]
    @test distsr[1] isa SampleDist

    ca3 = Makie.convert_arguments(SlabInterval, [Normal(0, 1), Normal(1, 1)])
    posv, distsv = ca3[1]
    @test posv == [1.0, 2.0]
    @test length(distsv) == 2
end
