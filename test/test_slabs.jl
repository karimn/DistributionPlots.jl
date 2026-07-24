using DistributionPlots
using DistributionPlots: slab_curve, asdist
using Distributions
using Test

@testset "slab_curve" begin
    a = asdist(Normal(0, 1))
    xs, th = slab_curve(a; kind=:pdf, n=201)
    @test length(xs) == length(th) == 201
    @test th[argmin(abs.(xs))] ≈ pdf(Normal(0,1), 0.0) atol=1e-2   # peak at 0
    @test all(≥(0), th)

    # cdf is monotone nondecreasing from ~0 to ~1
    xs2, cd = slab_curve(a; kind=:cdf, n=201)
    @test issorted(cd)
    @test cd[1] < 0.05 && cd[end] > 0.95

    # ccdf is 1 - cdf
    _, cc = slab_curve(a; kind=:ccdf, n=201)
    @test cc ≈ 1 .- cd atol=1e-6

    # histogram returns bin centers and counts/density
    s = asdist(randn(5000))
    xh, h = slab_curve(s; kind=:histogram)
    @test length(xh) == length(h)
    @test all(≥(0), h)
end
