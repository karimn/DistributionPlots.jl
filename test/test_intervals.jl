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

    # hdci on a symmetric sample ≈ qi (narrowest contiguous window)
    rh = point_interval(s; widths=[0.9], interval=:hdci)
    @test rh[1].interval == :hdci
    @test (rh[1].upper - rh[1].lower) ≤ (rmean[1].upper - rmean[1].lower) + 1e-9

    # Tables.jl compatibility: a Vector{NamedTuple} is a valid row table
    @test Tables.istable(rows)
    @test Tables.rowaccess(rows)
end
