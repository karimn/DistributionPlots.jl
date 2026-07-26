using DistributionPlots
using DistributionPlots: slab_polygon, interval_segment, point_marker, normalize_thickness,
                         interval_linewidths
# Makie is a dep of DistributionPlots but not a direct dep of test/Project.toml
# (only CairoMakie is); re-export its binding into Main here, same pattern used
# above for Tables/Statistics.
using DistributionPlots: Makie
using .Makie: Point2f
using Test

@testset "geometry" begin
    @test normalize_thickness([0.0, 1.0, 2.0], :all) == [0.0, 0.5, 1.0]
    @test normalize_thickness([0.0, 1.0, 2.0], :none) == [0.0, 1.0, 2.0]

    # :each normalizes to this vector's own max (the pre-fix :all behavior)
    @test normalize_thickness([0.0, 1.0, 2.0], :each) == [0.0, 0.5, 1.0]
    # :all with a supplied global divisor preserves RELATIVE height (peak need not be 1)
    @test normalize_thickness([0.0, 1.0, 2.0], :all; globalmax=4.0) == [0.0, 0.25, 0.5]
    # standalone :all still falls back to own max (backward compatible)
    @test normalize_thickness([0.0, 1.0, 2.0], :all) == [0.0, 0.5, 1.0]
    @test_throws ArgumentError normalize_thickness([1.0, 2.0], :bogus)

    xs = [0.0, 1.0, 2.0]
    th = [0.0, 1.0, 0.0]                        # normalised triangle
    # vertical, side=:top, position=5, scale=1 → polygon x ∈ [5, 6], y follows xs
    poly = slab_polygon(xs, th; position=5.0, orientation=:vertical,
                        side=:top, justification=0.0, scale=1.0, normalize=:none)
    @test poly isa Vector{Point2f}
    @test all(p -> p[1] ≥ 5.0 - 1e-6, poly)      # never left of baseline
    @test maximum(p -> p[1], poly) ≈ 6.0         # peak thickness 1 * scale 1
    @test any(p -> p[2] ≈ 1.0, poly)             # value axis is y

    # interval segment: vertical orientation → x fixed at position, y = [lo, hi]
    p1, p2 = interval_segment(-1.0, 2.0; position=5.0, orientation=:vertical)
    @test p1 == Point2f(5.0, -1.0)
    @test p2 == Point2f(5.0, 2.0)

    pm = point_marker(0.5; position=5.0, orientation=:vertical)
    @test pm == Point2f(5.0, 0.5)

    # horizontal orientation swaps axes
    p1h, p2h = interval_segment(-1.0, 2.0; position=5.0, orientation=:horizontal)
    @test p1h == Point2f(-1.0, 5.0)
    @test p2h == Point2f(2.0, 5.0)
end

@testset "interval_linewidths" begin
    # automatic: narrower width → thicker line, so credible levels nest visibly
    lws = interval_linewidths([0.66, 0.95], Makie.automatic)
    @test lws[0.66] > lws[0.95]

    # a single width still gets a sensible (non-zero) automatic weight
    lws1 = interval_linewidths([0.9], Makie.automatic)
    @test lws1[0.9] > 0

    # scalar behaves as before #5: every width gets the same weight
    lwss = interval_linewidths([0.66, 0.95], 4.0)
    @test lwss[0.66] == lwss[0.95] == 4.0

    # vector: one entry per distinct width, in widest→narrowest order
    lwsv = interval_linewidths([0.66, 0.95], [3.0, 9.0])
    @test lwsv[0.95] == 3.0   # widest gets the first entry
    @test lwsv[0.66] == 9.0  # narrowest gets the second entry

    # vector length must match the number of distinct widths
    @test_throws ArgumentError interval_linewidths([0.66, 0.95], [1.0, 2.0, 3.0])
end
