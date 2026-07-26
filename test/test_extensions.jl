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

@testset "RVars dimension names and labels" begin
    CairoMakie.activate!()
    # a[trial, arm]: arm 2 is shifted by +2 so pooling along each dimension is
    # checkable by the resulting medians rather than only by shape.
    nd, ntrial, narm = 400, 3, 2
    arr = randn(nd, ntrial, narm)
    arr[:, :, 2] .+= 2.0
    rv = RVar(arr; dimnames = (:trial, :arm), dimlabels = (nothing, ["control", "drug"]))

    @testset "rank ≥ 2 refuses to guess" begin
        # Before dimension support this silently took size(raw, 2) and vec'd a 2-D
        # slice, mixing the arms together into ntrial garbage distributions. It must
        # name the dimensions it found instead.
        err = try
            DistributionPlots._to_dist_args(rv); nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("2 dimensions", err.msg)
        @test occursin("trial", err.msg) && occursin("arm", err.msg)
    end

    @testset "dim= selects the axis dimension" begin
        posa, dista = DistributionPlots._to_dist_args(rv; dim = :arm)
        @test posa == [1.0, 2.0]
        @test length(dista) == 2
        # arm 1 pools to ~0, arm 2 to ~2 — the other axis is pooled, not interleaved
        meds = [point_interval(d)[1].value for d in dista]
        @test meds[1] ≈ 0.0 atol = 0.25
        @test meds[2] ≈ 2.0 atol = 0.25

        post, distt = DistributionPlots._to_dist_args(rv; dim = :trial)
        @test post == [1.0, 2.0, 3.0]
        @test length(distt) == 3
        # every trial pools both arms, so each is a two-component mixture near 1
        @test all(d -> isapprox(point_interval(d)[1].value, 1.0, atol = 0.5), distt)

        # a dimension may also be given by position
        @test DistributionPlots._to_dist_args(rv; dim = 2)[1] == posa
    end

    @testset "dim= is validated" begin
        @test_throws ArgumentError DistributionPlots._to_dist_args(rv; dim = :nope)
        @test_throws ArgumentError DistributionPlots._to_dist_args(rv; dim = 5)
        # naming a dimension of an RVar that has none should say so
        plain = RVar(randn(100, 2, 2))
        @test_throws ArgumentError DistributionPlots._to_dist_args(plain; dim = :arm)
        # a scalar RV has no dimension to select
        @test_throws ArgumentError DistributionPlots._to_dist_args(RVar(randn(100)); dim = :arm)
    end

    @testset "dimticks" begin
        @test dimticks(rv, :arm) == (1:2, ["control", "drug"])
        # an unlabelled axis falls back to its integer positions
        @test dimticks(rv, :trial) == (1:3, ["1", "2", "3"])
        @test dimticks(rv, 2) == (1:2, ["control", "drug"])
        # labels are stringified only at this boundary; non-String label types work
        rvi = RVar(randn(50, 3); dimnames = (:g,), dimlabels = ([10, 20, 30],))
        @test dimticks(rvi) == (1:3, ["10", "20", "30"])
        # element names of a named vector RV stand in when the axis is unlabelled
        rvn = RVar(randn(50, 2); names = [:a, :b])
        @test dimticks(rvn) == (1:2, ["a", "b"])
        @test_throws ArgumentError dimticks(RVar(randn(50)))
    end

    @testset "recipes accept dim=" begin
        for f in (halfeye, pointinterval, dotsinterval, interval, slab)
            @test f(rv; dim = :arm) isa Makie.FigureAxisPlot
        end
    end
end

using MCMCChains

@testset "MCMCChains extension" begin
    CairoMakie.activate!()
    chn = Chains(randn(500, 3, 2), [:a, :b, :c])       # 500 iters, 3 params, 2 chains
    @test (halfeye(chn) isa Makie.FigureAxisPlot)      # 3 positions, chains pooled
    pos, dists = DistributionPlots._to_dist_args(chn)
    @test pos == [1.0, 2.0, 3.0]

    # RVars ≥ 0.4 returns a NamedTuple of shaped RVars from RVar(chn); plotting a
    # whole fit is the flat "one distribution per scalar parameter" view, so the
    # extension must ask for flat=true rather than hand a NamedTuple to the converter.
    @test RVar(chn) isa NamedTuple
    @test length(dists) == 3
end

@testset "shaped parameters extracted from a Chains" begin
    CairoMakie.activate!()
    # a[trial, arm] arrives from a fit as flat a[1,1], a[2,1], ... columns; rvars
    # regroups them and attaches the names and labels the chain never recorded.
    names = [Symbol("a[$i,$j]") for j in 1:2 for i in 1:3]
    chn = Chains(randn(200, 6, 2), names)
    p = rvars(chn; dims = (a = (:trial, :arm),), labels = (arm = ["control", "drug"],))
    @test RVars.dimnames(p.a) == (:trial, :arm)
    @test RVars.dimlabels(p.a, :arm) == ["control", "drug"]
    @test dimticks(p.a, :arm) == (1:2, ["control", "drug"])
    @test halfeye(p.a; dim = :arm) isa Makie.FigureAxisPlot
end

@testset "two-argument (x, values) form" begin
    CairoMakie.activate!()
    x = repeat([1.0, 2.0, 3.0], inner = 400)
    y = randn(1200) .+ repeat([0.0, 2.0, 4.0], inner = 400)

    pos, dists = DistributionPlots._to_dist_args(x, y)
    @test pos == [1.0, 2.0, 3.0]
    @test length(dists) == 3
    # values are grouped by their position, not concatenated
    meds = [point_interval(d)[1].value for d in dists]
    @test meds ≈ [0.0, 2.0, 4.0] atol = 0.3

    # positions come back sorted and deduplicated regardless of input order
    shuffled = reverse(x)
    @test DistributionPlots._to_dist_args(shuffled, reverse(y))[1] == [1.0, 2.0, 3.0]

    @test_throws DimensionMismatch DistributionPlots._to_dist_args(x, y[1:10])

    for f in (halfeye, pointinterval, interval, dots, dotsinterval, slab)
        @test f(x, y) isa Makie.FigureAxisPlot
    end
end

@testset "dodge" begin
    CairoMakie.activate!()
    x = repeat([1.0, 2.0], inner = 300)
    y = randn(600)

    @test halfeye(x, y; dodge = 1, n_dodge = 2) isa Makie.FigureAxisPlot
    @test dotsinterval(x, y; dodge = 2, n_dodge = 2) isa Makie.FigureAxisPlot

    # AlgebraOfGraphics calls the recipe once per dodge group with a scalar `dodge`
    # and never sets `n_dodge`, so the group count cannot be inferred from one call.
    # Guessing would mis-size every group; demand it instead.
    @test_throws ArgumentError halfeye(x, y; dodge = 1)

    # placement: n groups share one slot, centred on it, and each narrows to its share
    s1, w1 = DistributionPlots.dodge_placement(1, 2, 0.0, 1.0)
    s2, w2 = DistributionPlots.dodge_placement(2, 2, 0.0, 1.0)
    @test s1 ≈ -0.25 && s2 ≈ 0.25
    @test w1 ≈ 0.5 && w2 ≈ 0.5
    @test s1 + s2 ≈ 0.0                      # symmetric about the position
    # a single group, or none at all, is left where it is at full width
    @test DistributionPlots.dodge_placement(1, 1, 0.03, 1.0) == (0.0, 1.0)
    @test DistributionPlots.dodge_placement(nothing, nothing, 0.03, 1.0) == (0.0, 1.0)
    # three groups each take a third of the slot, the middle one staying put
    s_mid, w_mid = DistributionPlots.dodge_placement(2, 3, 0.0, 1.0)
    @test s_mid ≈ 0.0 atol = 1e-12
    @test w_mid ≈ 1 / 3
    # gaps eat into the group width but keep the arrangement centred
    sa, wa = DistributionPlots.dodge_placement(1, 2, 0.1, 1.0)
    sb, _ = DistributionPlots.dodge_placement(2, 2, 0.1, 1.0)
    @test wa ≈ 0.45 && sa + sb ≈ 0.0
    # `width` scales the whole arrangement
    @test DistributionPlots.dodge_placement(1, 2, 0.0, 0.5)[1] ≈ -0.125

    @test_throws ArgumentError DistributionPlots.dodge_placement(3, 2, 0.0, 1.0)
    # one plot call draws one group, so a per-element dodge vector is rejected
    @test_throws ArgumentError DistributionPlots.dodge_placement([1, 2], 2, 0.0, 1.0)
end

using AlgebraOfGraphics
const AoG = AlgebraOfGraphics

@testset "AlgebraOfGraphics extension" begin
    CairoMakie.activate!()
    # Single mapped column: our one-argument convert_arguments path, equivalent to
    # calling interval(randn(1500)) outside AoG.
    spec = AoG.data((y = randn(1500),)) * AoG.mapping(:y) * AoG.visual(Interval)
    fg = AoG.draw(spec)
    # `FigureGrid` is defined by AlgebraOfGraphics itself (not re-exported through
    # Makie), so it's asserted qualified as `AoG.FigureGrid` here.
    @test fg isa AoG.FigureGrid
end

@testset "AlgebraOfGraphics: labels drive colour, facets and legend" begin
    CairoMakie.activate!()
    nd, ntrial, narm = 200, 3, 2
    arr = randn(nd, ntrial, narm)
    arr[:, :, 2] .+= 2.0
    rv = RVar(arr; nchains = 2, dimnames = (:trial, :arm),
              dimlabels = (nothing, ["control", "drug"]))

    # RVars.gather_draws is the bridge: a labelled array becomes a long table whose
    # columns are named after the dimensions and hold the categories themselves.
    tbl = RVars.gather_draws(rv)
    @test keys(tbl) == (:trial, :arm, :chain, :draw, :value)
    @test length(tbl.value) == nd * ntrial * narm
    @test sort(unique(tbl.arm)) == ["control", "drug"]

    nlegends(fg) = count(x -> x isa Makie.Legend, fg.figure.content)
    naxes(fg) = count(x -> x isa Makie.Axis, fg.figure.content)

    @testset "colour + facet" begin
        fg = AoG.draw(AoG.data(tbl) *
                      AoG.mapping(:arm, :value; color = :arm, layout = :trial) *
                      AoG.visual(HalfEye))
        @test fg isa AoG.FigureGrid
        @test naxes(fg) == ntrial              # one panel per trial
        @test nlegends(fg) == 1                # legend_elements makes this possible
        ax = fg.figure.content[1]
        # AoG maps the categorical column to codes and sets the ticks itself, so the
        # axis speaks in the model's own terms without any help from the recipe.
        @test ax.xticks[] == (Base.OneTo(2), ["control", "drug"])
    end

    @testset "slab_color is separately mappable (ggdist's fill)" begin
        fg = AoG.draw(AoG.data(tbl) *
                      AoG.mapping(:arm, :value; slab_color = :arm) *
                      AoG.visual(HalfEye))
        @test fg isa AoG.FigureGrid
        @test nlegends(fg) == 1
    end

    @testset "dodge" begin
        fg = AoG.draw(AoG.data(tbl) *
                      AoG.mapping(:trial, :value; color = :arm, dodge = :arm) *
                      AoG.visual(HalfEye; n_dodge = 2))
        @test fg isa AoG.FigureGrid
        @test nlegends(fg) == 1
    end

    @testset "every public recipe type registers its aesthetics and legend" begin
        for T in (SlabInterval, HalfEye, Eye, CcdfInterval, CdfInterval,
                  GradientInterval, HistInterval, Slab, Interval, PointInterval,
                  Spike, Dots, DotsInterval)
            fg = AoG.draw(AoG.data(tbl) *
                          AoG.mapping(:arm, :value; color = :arm) * AoG.visual(T))
            @test fg isa AoG.FigureGrid
            @test nlegends(fg) == 1
        end
    end
end
