using KernelDensity: kde, pdf as kdepdf

"""
    kde_reflected(samples; bounds, npoints=512)

Kernel density estimate that reflects samples across any *finite* bound so the
density is not underestimated near the boundary (ggdist's `density = "bounded"`).
Returns `(xs, dens)` evaluated on `npoints` points spanning `bounds`.
"""
function kde_reflected(samples::AbstractVector{<:Real};
                       bounds::Tuple{<:Real,<:Real},
                       npoints::Int=512)
    lo, hi = float(bounds[1]), float(bounds[2])
    data = collect(Float64, samples)

    reflect_lo = isfinite(lo)
    reflect_hi = isfinite(hi)
    augmented = copy(data)
    reflect_lo && append!(augmented, 2lo .- data)   # mirror across lower bound
    reflect_hi && append!(augmented, 2hi .- data)   # mirror across upper bound

    k = kde(augmented)                                # KernelDensity picks bandwidth
    xs = collect(range(lo, hi; length=npoints))
    dens = kdepdf(k, xs)
    # Reflection triples (or doubles) total mass inside [lo,hi]; renormalise so
    # the curve integrates to ~1 over the true support.
    scale = 1 + reflect_lo + reflect_hi
    dens .*= scale
    dens .= max.(dens, 0.0)
    return xs, dens
end

density_at(d::AnalyticDist, xs::AbstractVector{<:Real}; bounds=support(d)) =
    pdf.(d.dist, xs)

function density_at(d::SampleDist, xs::AbstractVector{<:Real}; bounds=support(d))
    gx, gd = kde_reflected(d.samples; bounds=bounds, npoints=max(512, length(xs)))
    # linear-interpolate the grid density onto the requested xs
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
