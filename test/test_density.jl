using DistributionPlots
using DistributionPlots: asdist, density_at, kde_reflected, support, slab_curve
using Distributions
using Test

@testset "density" begin
    # Analytic density_at is exact pdf
    a = asdist(Normal(0, 1))
    @test density_at(a, [0.0]) ≈ [pdf(Normal(0,1), 0.0)]

    # Reflected KDE on [0,1]-bounded samples: mass near a boundary is not
    # underestimated the way an unbounded KDE would underestimate it.
    Random_seeded = Beta(2, 8)
    samples = quantile.(Beta(2,8), range(0.001, 0.999; length=4000))  # deterministic
    xs, dens = kde_reflected(samples; bounds=(0.0, 1.0), npoints=512)
    @test all(≥(0), dens)
    @test minimum(xs) ≥ -1e-9 && maximum(xs) ≤ 1 + 1e-9
    # integral ≈ 1 over [0,1]
    area = sum((dens[2:end] .+ dens[1:end-1]) ./ 2 .* diff(xs))
    @test area ≈ 1.0 atol=0.05

    # Reflected density near the lower boundary exceeds the naive-unbounded
    # estimate at the same point (reflection adds the mirrored mass).
    s = asdist(samples)
    d_lo = density_at(s, [0.05]; bounds=(0.0, 1.0))[1]
    @test d_lo > 0

    # Unbounded samples (the default): the slab window extends PAST the data
    # extremes so the KDE tapers naturally — it is NOT reflected/cut at min/max.
    usamp = quantile.(Normal(0, 1), range(0.0005, 0.9995; length=4000))
    uxs, uth = slab_curve(usamp; kind=:pdf, n=201)
    @test minimum(uxs) < minimum(usamp) - 1e-9      # window padded below the data
    @test maximum(uxs) > maximum(usamp) + 1e-9      # window padded above the data
    @test uth[1]   < 0.05 * maximum(uth)            # far-tail density tapers to ~0
    @test uth[end] < 0.05 * maximum(uth)
    @test uth[argmin(abs.(uxs))] ≈ pdf(Normal(0,1), 0.0) atol=0.03   # peak ≈ N(0,1) pdf
end
