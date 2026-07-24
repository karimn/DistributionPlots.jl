using DistributionPlots
using DistributionPlots: dot_layout, asdist
using Distributions
using Test

@testset "dot_layout" begin
    a = asdist(Normal(0, 1))
    lay = dot_layout(a; ndots=50)
    @test length(lay.x) == 50
    @test length(lay.y) == 50
    @test all(≥(1), lay.y)                       # stacks are 1-based
    @test lay.binwidth > 0
    # dots are sorted along x (quantile dotplot)
    @test issorted(lay.x)
    # near the mode (x≈0) stacks are taller than in the tails
    center_heights = maximum(lay.y[abs.(lay.x) .< 0.5])
    tail_heights = maximum(lay.y[abs.(lay.x) .> 2.0])
    @test center_heights ≥ tail_heights
end
