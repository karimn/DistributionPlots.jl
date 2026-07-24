using DistributionPlots
using DistributionPlots: point_interval, point_summary, asdist
using Test
# mean, median, quantile, and Tables are brought into Main by test/runtests.jl

@testset "intervals: qi + hdci" begin
    s = collect(1.0:100.0)              # uniform-ish sample
    rows = point_interval(s; widths=[0.5, 0.9], point=:median, interval=:qi)
    @test length(rows) == 2
    @test rows[1].width == 0.5
    @test rows[1].point == :median
    @test rows[1].interval == :qi
    @test rows[1].value ≈ median(s)
    # 50% qi = [q25, q75]
    @test rows[1].lower ≈ quantile(s, 0.25)
    @test rows[1].upper ≈ quantile(s, 0.75)
    # 90% qi = [q05, q95]
    @test rows[2].lower ≈ quantile(s, 0.05)
    @test rows[2].upper ≈ quantile(s, 0.95)

    # mean point summary
    rmean = point_interval(s; widths=[0.9], point=:mean)
    @test rmean[1].value ≈ mean(s)

    # hdci on a symmetric sample ≈ qi (narrowest contiguous window).
    # `_hdci` uses a type-5 empirical quantile function (matching ggdist's
    # `hdci`) while `:qi` uses a type-7 quantile (matching ggdist's `qi`), so
    # for this perfectly evenly-spaced integer grid the two conventions'
    # boundary interpolation can disagree by up to one grid step (here
    # 89.5 vs 89.1); allow that margin rather than requiring bit-exactness
    # across two deliberately-different quantile conventions.
    rh = point_interval(s; widths=[0.9], interval=:hdci)
    @test rh[1].interval == :hdci
    @test (rh[1].upper - rh[1].lower) ≤ (rmean[1].upper - rmean[1].lower) + 1.0

    # Tables.jl compatibility: a Vector{NamedTuple} is a valid row table
    @test Tables.istable(rows)
    @test Tables.rowaccess(rows)
end

@testset "intervals: mode + disjoint hdi" begin
    using Distributions
    # Bimodal mixture: two well-separated Normals → hdi returns TWO pieces at 0.95.
    mix = MixtureModel([Normal(-3, 0.4), Normal(3, 0.4)], [0.5, 0.5])
    draws = reduce(vcat, (rand(Normal(-3,0.4), 5000), rand(Normal(3,0.4), 5000)))
    rows = point_interval(draws; widths=[0.95], interval=:hdi, point=:median)
    @test length(rows) ≥ 2                       # disjoint → multiple rows
    @test all(r -> r.interval == :hdi && r.width == 0.95, rows)
    # the two pieces straddle the two modes
    los = sort(getproperty.(rows, :lower))
    @test los[1] < 0 && los[end] > 0

    # mode of a right-skewed sample is left of the mean
    skew = rand(Exponential(1.0), 20000)
    rmode = point_interval(skew; widths=[0.9], point=:mode)
    @test rmode[1].value < mean(skew)
end
