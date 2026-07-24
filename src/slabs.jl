using StatsBase: fit, Histogram

"""
    slab_curve(d; kind=:pdf, n=201, trim=0.001)

Return `(xs, thickness)` describing the slab shape. `kind`:
- `:pdf`  — density (KDE for samples, exact pdf for analytic)
- `:cdf`  — cumulative distribution
- `:ccdf` — complementary cdf (`1 - cdf`)
- `:histogram` — binned counts (samples only)
"""
function slab_curve(d; kind::Symbol=:pdf, n::Int=201, trim::Real=0.001)
    dist = asdist(d)
    lo, hi = support(dist; trim=trim)
    if kind === :histogram
        return _histogram_curve(dist)
    end
    xs = collect(range(lo, hi; length=n))
    if kind === :pdf
        return xs, density_at(dist, xs; bounds=(lo, hi))
    elseif kind === :cdf
        return xs, [cdf_at(dist, x) for x in xs]
    elseif kind === :ccdf
        return xs, [1 - cdf_at(dist, x) for x in xs]
    else
        throw(ArgumentError("slab_curve: unknown kind :$kind (expected :pdf, :cdf, :ccdf, or :histogram)"))
    end
end

function _histogram_curve(d::SampleDist)
    h = fit(Histogram, d.samples; nbins=_fd_nbins(d.samples))
    edges = collect(h.edges[1])
    centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
    weights = float.(h.weights)
    return centers, weights
end

_histogram_curve(::AnalyticDist) =
    throw(ArgumentError("slab_curve: kind=:histogram requires samples, not an analytic distribution"))

# Freedman–Diaconis bin count
function _fd_nbins(x::Vector{Float64})
    n = length(x)
    iqr = quantile(x, 0.75) - quantile(x, 0.25)
    bw = iqr > 0 ? 2 * iqr / cbrt(n) : (maximum(x) - minimum(x)) / max(1, round(Int, sqrt(n)))
    bw ≤ 0 && return 1
    return max(1, ceil(Int, (maximum(x) - minimum(x)) / bw))
end
