using DistributionPlots
using DistributionPlots: SampleDist, AnalyticDist, asdist, support, quantile_at, cdf_at
using Distributions
using Test

@testset "interface" begin
    s = asdist([1.0, 2.0, 3.0, 4.0])
    @test s isa SampleDist
    @test quantile_at(s, 0.5) ≈ 2.5          # type-7 median of 1:4
    @test cdf_at(s, 2.5) ≈ 0.5
    @test support(s) == (1.0, 4.0)

    a = asdist(Normal(0, 1))
    @test a isa AnalyticDist
    @test quantile_at(a, 0.5) ≈ 0.0 atol=1e-12
    @test cdf_at(a, 0.0) ≈ 0.5
    lo, hi = support(a; trim=0.001)
    @test lo ≈ quantile(Normal(0,1), 0.001)
    @test hi ≈ quantile(Normal(0,1), 0.999)

    b = asdist(Uniform(0, 1))                # finite support: not trimmed
    @test support(b) == (0.0, 1.0)

    @test_throws ArgumentError asdist(Float64[])            # empty rejected
    @test_logs (:warn,) asdist([1.0, NaN, 3.0])            # NaN dropped with warn
    @test asdist([1.0, NaN, 3.0]).samples == [1.0, 3.0]
end
