using DistributionPlots
using DistributionPlots: SlabInterval, _to_dist_args, AnalyticDist, SampleDist
using DistributionPlots: Makie
using .Makie
using Distributions
using Test
using CairoMakie

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

@testset "slabinterval renders" begin
    CairoMakie.activate!()
    f = slabinterval(randn(2000))
    @test f isa Makie.FigureAxisPlot
    ax = f.axis
    # at least one poly (slab) and some line/scatter children exist
    @test !isempty(ax.scene.plots)

    # toggles: show only the interval
    f2 = slabinterval(randn(2000); show_slab=false, show_point=false)
    @test f2 isa Makie.FigureAxisPlot

    # saving to a headless PNG does not error
    tmp = tempname() * ".png"
    save(tmp, f)
    @test isfile(tmp)

    # multi-width nesting (#5): the default is automatic per-width line weight,
    # not one flat linewidth for every credible level
    th = Makie.default_theme(nothing, SlabInterval)
    @test th[:interval_linewidth][] === Makie.automatic
    @test th[:point_size][] > 10   # bumped up for contrast against the interval

    # a scalar still means "every width the same", pre-#5 behaviour
    @test (slabinterval(randn(1500); interval_linewidth=4.0) isa Makie.FigureAxisPlot)
    # a vector (one entry per distinct width) draws without error
    @test (slabinterval(randn(1500); widths=[0.66, 0.95],
                        interval_linewidth=[3.0, 9.0]) isa Makie.FigureAxisPlot)
    # wrong-length vector is rejected loudly rather than silently mismatched
    @test_throws ArgumentError slabinterval(randn(1500); widths=[0.66, 0.95],
                                            interval_linewidth=[1.0, 2.0, 3.0])
end

@testset "children + pre-summarised" begin
    CairoMakie.activate!()
    for f in (halfeye, eye, ccdfinterval, cdfinterval, gradientinterval,
              histinterval, slab, interval, pointinterval, spike)
        @test (f(randn(1500)) isa Makie.FigureAxisPlot)
    end

    # halfeye defaults: slab on top, point shown
    th = Makie.default_theme(nothing, HalfEye)
    @test th[:side][] == :top

    # override a child default
    @test (halfeye(randn(1500); interval=:hdi) isa Makie.FigureAxisPlot)

    # pre-summarised pointinterval: positions + point + lower/upper
    @test (pointinterval([1.0,2.0], [0.0,1.0], [-1.0,0.0], [1.0,2.0]) isa Makie.FigureAxisPlot)

    # slab from pre-summarised data is rejected loudly (not a bare MethodError)
    @test_throws ArgumentError DistributionPlots._reject_slab_summary()
end

@testset "dots" begin
    CairoMakie.activate!()
    @test (dots(randn(1000)) isa Makie.FigureAxisPlot)
    @test (dotsinterval(randn(1000)) isa Makie.FigureAxisPlot)
end

@testset "orientation" begin
    CairoMakie.activate!()

    # default orientation is :vertical for every family that carries the attribute
    for R in (SlabInterval, Dots, DotsInterval)
        th = Makie.default_theme(nothing, R)
        @test haskey(th, :orientation)
        @test th[:orientation][] == :vertical
    end

    # an unrecognised orientation throws naming the valid options rather than
    # silently falling back to :vertical
    @test_throws ArgumentError halfeye(randn(500); orientation=:sideways)
    @test_throws ArgumentError dots(randn(500); orientation=:sideways)

    # renders without error under :horizontal, for both families
    @test (halfeye(randn(1500); orientation=:horizontal) isa Makie.FigureAxisPlot)
    @test (pointinterval(randn(1500); orientation=:horizontal) isa Makie.FigureAxisPlot)
    @test (dots(randn(500); orientation=:horizontal) isa Makie.FigureAxisPlot)
    @test (dotsinterval(randn(500); orientation=:horizontal) isa Makie.FigureAxisPlot)

    # a horizontal plot is the vertical one with coordinates swapped — checked on
    # actual point coordinates via the pre-summarised `pointinterval`, whose
    # linesegments!/scatter! calls land directly on the returned axis (not nested
    # inside a recipe), matching the pattern used elsewhere in this file.
    _, axv, _ = pointinterval([1.0, 2.0], [0.0, 1.0], [-1.0, 0.0], [1.0, 2.0])
    _, axh, _ = pointinterval([1.0, 2.0], [0.0, 1.0], [-1.0, 0.0], [1.0, 2.0]; orientation=:horizontal)

    lines_v = filter(p -> p isa Makie.LineSegments, axv.scene.plots)
    lines_h = filter(p -> p isa Makie.LineSegments, axh.scene.plots)
    @test length(lines_v) == length(lines_h) == 2
    for (lv, lh) in zip(lines_v, lines_h)
        segs_v, segs_h = lv[1][], lh[1][]
        @test [Makie.Point2f(pt[2], pt[1]) for pt in segs_v] == segs_h
    end

    scats_v = filter(p -> p isa Makie.Scatter, axv.scene.plots)
    scats_h = filter(p -> p isa Makie.Scatter, axh.scene.plots)
    @test length(scats_v) == length(scats_h) == 2
    for (sv, sh) in zip(scats_v, scats_h)
        pv, ph = sv[1][][1], sh[1][][1]
        @test Makie.Point2f(pv[2], pv[1]) == ph
    end
end

@testset "lineribbon" begin
    CairoMakie.activate!()
    xgrid = collect(0.0:0.5:10.0)
    # a growing-uncertainty fan: samples per x
    perx = [randn(800) .* (0.2 + 0.1x) .+ sin(x) for x in xgrid]
    @test (lineribbon(xgrid, perx) isa Makie.FigureAxisPlot)
end
