using KernelDensity: kde, pdf as kdepdf
using Statistics: std

# Silverman rule-of-thumb bandwidth; used to pad the slab window past the data
# so an unbounded KDE tapers to ~0 instead of being cut at the sample range.
function _kde_bandwidth(samples::AbstractVector{<:Real})
    x = collect(Float64, samples)
    n = length(x)
    s = n > 1 ? std(x) : 0.0
    return s > 0 ? 1.06 * s * n^(-1/5) : 1.0
end

"""
    kde_reflected(samples; bounds, eval=bounds, npoints=512)

Kernel density estimate that reflects samples across any *finite* `bounds` so the
density is not underestimated near a known boundary (ggdist's `density = "bounded"`).
`bounds` are the reflection bounds — an infinite side is not reflected. The curve is
evaluated on `npoints` points spanning `eval` (defaults to `bounds`); density beyond a
finite reflection bound is set to zero. Returns `(xs, dens)`.
"""
function kde_reflected(samples::AbstractVector{<:Real};
                       bounds::Tuple{<:Real,<:Real},
                       eval::Tuple{<:Real,<:Real}=bounds,
                       npoints::Int=512)
    lo, hi = float(bounds[1]), float(bounds[2])
    data = collect(Float64, samples)

    reflect_lo = isfinite(lo)
    reflect_hi = isfinite(hi)
    augmented = copy(data)
    reflect_lo && append!(augmented, 2lo .- data)   # mirror across lower bound
    reflect_hi && append!(augmented, 2hi .- data)   # mirror across upper bound

    k = kde(augmented)                                # KernelDensity picks bandwidth
    xs = collect(range(float(eval[1]), float(eval[2]); length=npoints))
    dens = kdepdf(k, xs)
    # Reflection multiplies total mass inside the true support; renormalise.
    scale = 1 + reflect_lo + reflect_hi
    dens .*= scale
    # Zero out anything beyond a finite reflection bound (outside true support).
    @inbounds for i in eachindex(xs)
        (xs[i] < lo || xs[i] > hi) && (dens[i] = 0.0)
    end
    dens .= max.(dens, 0.0)
    return xs, dens
end

density_at(d::AnalyticDist, xs::AbstractVector{<:Real}; bounds=nothing) =
    pdf.(d.dist, xs)

function density_at(d::SampleDist, xs::AbstractVector{<:Real}; bounds=d.bounds)
    e = extrema(xs)
    gx, gd = kde_reflected(d.samples; bounds=bounds, eval=e, npoints=max(512, length(xs)))
    return _interp(gx, gd, xs)
end

function _interp(gx::Vector{Float64}, gd::Vector{Float64}, xs)
    out = similar(collect(Float64, xs))
    for (i, x) in enumerate(xs)
        if x ≤ gx[1]
            out[i] = gd[1]
        elseif x ≥ gx[end]
            out[i] = gd[end]
        else
            j = searchsortedlast(gx, x)
            t = (x - gx[j]) / (gx[j+1] - gx[j])
            out[i] = (1 - t) * gd[j] + t * gd[j+1]
        end
    end
    return out
end
