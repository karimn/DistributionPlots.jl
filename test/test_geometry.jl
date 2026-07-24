using DistributionPlots
using DistributionPlots: slab_polygon, interval_segment, point_marker, normalize_thickness
# Makie is a dep of DistributionPlots but not a direct dep of test/Project.toml
# (only CairoMakie is); re-export its binding into Main here, same pattern used
# above for Tables/Statistics.
using DistributionPlots: Makie
using .Makie: Point2f
using Test

@testset "geometry" begin
    @test normalize_thickness([0.0, 1.0, 2.0], :all) == [0.0, 0.5, 1.0]
    @test normalize_thickness([0.0, 1.0, 2.0], :none) == [0.0, 1.0, 2.0]

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
