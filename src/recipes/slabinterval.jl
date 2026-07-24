using Makie

# Shared converter: normalise any input to (positions, Vector{AbstractDist}).
# A single distribution/sample-vector is placed at position 1. Extensions
# (RVars, MCMCChains) add methods that expand multi-element inputs to 1:k.
_to_dist_args(x::Distributions.UnivariateDistribution) = ([1.0], AbstractDist[asdist(x)])
_to_dist_args(x::AbstractVector{<:Real}) = ([1.0], AbstractDist[asdist(x)])
function _to_dist_args(xs::AbstractVector{<:Distributions.UnivariateDistribution})
    return (collect(1.0:length(xs)), AbstractDist[asdist(x) for x in xs])
end

@recipe(SlabInterval, data) do scene
    Attributes(
        slab_type = :pdf,          # :pdf, :cdf, :ccdf, :histogram
        interval = :qi,            # :qi, :hdci, :hdi
        point = :median,           # :mean, :median, :mode
        widths = [0.66, 0.95],
        side = :top,               # :top, :bottom, :both
        justification = 0.0,
        scale = 0.9,
        normalize = :all,          # :all, :none
        show_slab = true,
        show_interval = true,
        show_point = true,
        trim = 0.001,
        n = 201,
        color = :black,
        colormap = :blues,
        colorrange = Makie.automatic,
        interval_linewidth = 6,
        point_size = 10,
    )
end

# Register convert_arguments for the base recipe. Task 12 loops this over all
# public recipe types; the extensions add per-type methods for RVar/Chains.
Makie.convert_arguments(::Type{<:SlabInterval}, x::Distributions.UnivariateDistribution) =
    (_to_dist_args(x),)
Makie.convert_arguments(::Type{<:SlabInterval}, x::AbstractVector{<:Real}) =
    (_to_dist_args(x),)
Makie.convert_arguments(::Type{<:SlabInterval}, xs::AbstractVector{<:Distributions.UnivariateDistribution}) =
    (_to_dist_args(xs),)

function Makie.plot!(p::SlabInterval)
    # p.data already resolves to the (positions, dists) tuple produced by
    # convert_arguments (Makie unwraps the recipe's single positional arg).
    args = p.data

    lift(args, p.slab_type, p.interval, p.point, p.widths, p.side, p.justification,
         p.scale, p.normalize, p.trim, p.n, p.show_slab, p.show_interval, p.show_point,
         p.color, p.interval_linewidth, p.point_size) do (positions, dists), st, iv, pt,
            ws, side, just, sc, nrm, tr, n, sslab, sint, spoint, col, ilw, psz

        for (pos, d) in zip(positions, dists)
            if sslab
                xs, th = slab_curve(d; kind=st, n=n, trim=tr)
                poly!(p, slab_polygon(xs, th; position=pos, orientation=:vertical,
                                      side=side, justification=just, scale=sc, normalize=nrm);
                      color=col)
            end
            rows = point_interval(d; widths=ws, point=pt, interval=iv)
            if sint
                for r in rows
                    a, b = interval_segment(r.lower, r.upper; position=pos, orientation=:vertical)
                    linesegments!(p, [a, b]; linewidth=ilw, color=col)
                end
            end
            if spoint && !isempty(rows)
                scatter!(p, [point_marker(rows[1].value; position=pos, orientation=:vertical)];
                         markersize=psz, color=col)
            end
        end
    end
    return p
end
