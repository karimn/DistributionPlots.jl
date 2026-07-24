abstract type AbstractDist end

struct SampleDist <: AbstractDist
    samples::Vector{Float64}
end

struct AnalyticDist{D<:Distributions.UnivariateDistribution} <: AbstractDist
    dist::D
end

"""
    asdist(x)

Normalise an input into an `AbstractDist`. Vectors of samples become
`SampleDist` (NaN/missing dropped with a warning; empty is an error);
`UnivariateDistribution`s become `AnalyticDist`.
"""
function asdist(x::AbstractVector{<:Real})
    clean = collect(Float64, Iterators.filter(!isnan, skipmissing(x)))
    ndropped = length(x) - length(clean)
    ndropped > 0 && @warn "asdist: dropped $ndropped NaN/missing value(s) from samples"
    isempty(clean) && throw(ArgumentError("asdist: no finite samples (got $(length(x)) values, all NaN/missing/empty)"))
    return SampleDist(clean)
end

asdist(d::Distributions.UnivariateDistribution) = AnalyticDist(d)
asdist(d::AbstractDist) = d

quantile_at(d::SampleDist, p::Real) = quantile(d.samples, p)          # StatsBase type-7
quantile_at(d::AnalyticDist, p::Real) = quantile(d.dist, p)

cdf_at(d::SampleDist, x::Real) = count(≤(x), d.samples) / length(d.samples)
cdf_at(d::AnalyticDist, x::Real) = cdf(d.dist, x)

function support(d::SampleDist; trim::Real=0.001)
    return (minimum(d.samples), maximum(d.samples))
end

function support(d::AnalyticDist; trim::Real=0.001)
    lo = minimum(d.dist)
    hi = maximum(d.dist)
    isfinite(lo) || (lo = quantile(d.dist, trim))
    isfinite(hi) || (hi = quantile(d.dist, 1 - trim))
    return (float(lo), float(hi))
end
