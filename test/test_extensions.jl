using DistributionPlots
using DistributionPlots: Makie
using .Makie
using RVars
using CairoMakie
using Test

@testset "RVars extension" begin
    CairoMakie.activate!()
    # scalar RV: 1000 draws → one distribution at position 1
    rv0 = RVar(randn(1000))
    pos, dists = DistributionPlots._to_dist_args(rv0)
    @test pos == [1.0]
    @test length(dists) == 1

    # vector RV: 3 elements → positions 1:3
    rv1 = RVar(randn(1000, 3))
    pos3, dists3 = DistributionPlots._to_dist_args(rv1)
    @test pos3 == [1.0, 2.0, 3.0]
    @test length(dists3) == 3

    # recipes accept a RVar directly
    @test (halfeye(rv1) isa Makie.FigureAxisPlot)
    @test (pointinterval(rv1) isa Makie.FigureAxisPlot)
end

using MCMCChains

@testset "MCMCChains extension" begin
    CairoMakie.activate!()
    chn = Chains(randn(500, 3, 2), [:a, :b, :c])       # 500 iters, 3 params, 2 chains
    @test (halfeye(chn) isa Makie.FigureAxisPlot)      # 3 positions, chains pooled
    pos, dists = DistributionPlots._to_dist_args(chn)
    @test pos == [1.0, 2.0, 3.0]
end

using AlgebraOfGraphics
const AoG = AlgebraOfGraphics

@testset "AlgebraOfGraphics extension" begin
    CairoMakie.activate!()
    # visual() over our recipe type composes and draws without error.
    #
    # Scoped to a single mapped column rather than the two-column (x, y)
    # form: our recipes' `convert_arguments` contract takes one positional
    # argument (a distribution, a sample vector, or a vector of
    # distributions/positions+dists for the pre-summarised child recipes),
    # not a pair of same-length raw (x, y) sample columns the way a generic
    # Makie `Scatter`-like recipe does. AlgebraOfGraphics (v0.13.1, per
    # test/Manifest.toml) is happy to call `Interval(y_column)` under
    # `mapping(:y)`, which round-trips through our single-arg
    # `convert_arguments(x::AbstractVector{<:Real})` path exactly as
    # `interval(randn(1500))` would outside of AoG.
    spec = AoG.data((y = randn(1500),)) *
           AoG.mapping(:y) *
           AoG.visual(Interval)
    fg = AoG.draw(spec)
    # `FigureGrid` is defined by AlgebraOfGraphics itself (not re-exported
    # through Makie), so it's asserted qualified as `AoG.FigureGrid` here.
    @test fg isa AoG.FigureGrid
end
