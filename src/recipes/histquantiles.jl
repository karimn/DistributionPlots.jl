using Makie

# nested central bands (lower, upper, alpha), widest and faintest first
const HISTQUANTILE_BANDS = ((0.1, 0.9, 0.2), (0.2, 0.8, 0.35), (0.3, 0.7, 0.5), (0.4, 0.6, 0.7))

"""
    histquantiles(vals; nbins = 40, lims = nothing, reference = nothing)

Betancourt-style histogram quantiles. `vals` is a draws × units matrix: each draw's
units are binned into a density histogram on shared bin edges, then the 10–90, 20–80,
30–70 and 40–60 % quantile bands across draws (and the median) are drawn per bin,
stepped to follow the bars.

`lims = (lo, hi)` fixes the binned range; by default it is the 0.2–99.8 % quantile
range of `vals`. Values outside the range are dropped, but the density is still
normalised by all units, so it integrates to the share of units inside the range.

`reference = (x, M)` overlays a per-draw reference density — `M` is draws × `length(x)`,
evaluated at `x` — as a 10–90 % band plus median. For a transformed axis, transform
`vals` and `lims` beforehand and give `reference` as a density on that same scale
(including any Jacobian).
"""
@recipe(HistQuantiles, vals) do scene
    Attributes(
        nbins = 40,
        lims = nothing,
        reference = nothing,
        color = :firebrick,
        mediancolor = :darkred,
        medianwidth = 1.5,
        referencecolor = :steelblue,
        referencewidth = 2,
    )
end

"Bin edges and a draws × bins density matrix of `vals` (draws × units)."
function hist_density(vals::AbstractMatrix{<:Real}, nbins::Integer, lims)
    lo, hi = isnothing(lims) ? quantile(filter(isfinite, vec(vals)), (0.002, 0.998)) : lims
    edges = collect(range(Float64(lo), Float64(hi); length = nbins + 1))
    w = (edges[end] - edges[1]) / nbins
    D, P = size(vals)
    H = zeros(D, nbins)
    for d in 1:D, v in view(vals, d, :)
        b = floor(Int, (v - edges[1]) / w) + 1
        v == edges[end] && (b = nbins)            # right edge is inclusive
        1 <= b <= nbins && (H[d, b] += 1)
    end
    H ./= P * w
    return edges, H
end

"Quantile `p` of each column of `M`."
colquantile(M, p) = [quantile(view(M, :, k), p) for k in axes(M, 2)]

"Stepped x/y for a per-bin quantity, so bands follow the histogram bars."
stepped(edges, y) = (repeat(edges; inner = 2)[2:end-1], repeat(y; inner = 2))

function Makie.plot!(p::HistQuantiles)
    hist = lift(p.vals, p.nbins, p.lims) do vals, nbins, lims
        hist_density(vals, nbins, lims)
    end
    for (k, (a, b, α)) in enumerate(HISTQUANTILE_BANDS)
        pts = lift(hist) do (edges, H)
            x, lo = stepped(edges, colquantile(H, a))
            _, hi = stepped(edges, colquantile(H, b))
            (x, lo, hi)
        end
        band!(p, lift(first, pts), lift(t -> t[2], pts), lift(t -> t[3], pts);
              color = lift(c -> (c, α), p.color))
    end
    med = lift(h -> stepped(h[1], colquantile(h[2], 0.5)), hist)
    lines!(p, lift(first, med), lift(last, med); color = p.mediancolor, linewidth = p.medianwidth)

    # empty vectors when there is no reference keep the plot objects (and layout) stable
    ref = lift(p.reference) do r
        isnothing(r) && return (Float64[], Float64[], Float64[], Float64[])
        x, M = r
        size(M, 2) == length(x) ||
            throw(ArgumentError("reference matrix must have one column per x ($(length(x))), got $(size(M, 2))"))
        (collect(Float64, x), colquantile(M, 0.1), colquantile(M, 0.9), colquantile(M, 0.5))
    end
    band!(p, lift(r -> r[1], ref), lift(r -> r[2], ref), lift(r -> r[3], ref);
          color = lift(c -> (c, 0.25), p.referencecolor))
    lines!(p, lift(r -> r[1], ref), lift(r -> r[4], ref);
           color = p.referencecolor, linewidth = p.referencewidth)
    return p
end

Makie.convert_arguments(::Type{<:HistQuantiles}, vals::AbstractMatrix{<:Real}) = (Matrix{Float64}(vals),)
