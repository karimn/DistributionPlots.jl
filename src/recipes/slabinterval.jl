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
        # ggdist parity: the slab is a grey, semi-transparent fill so the
        # black point/interval drawn on top of it stays visible instead of
        # being occluded (both would be solid black otherwise).
        slab_color = :grey55,
        slab_alpha = 0.7,
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
         p.color, p.slab_color, p.slab_alpha, p.interval_linewidth, p.point_size) do (positions, dists), st, iv, pt,
            ws, side, just, sc, nrm, tr, n, sslab, sint, spoint, col, scol, salpha, ilw, psz

        for (pos, d) in zip(positions, dists)
            # slab drawn first so the point/interval below render on top of it
            if sslab
                xs, th = slab_curve(d; kind=st, n=n, trim=tr)
                poly!(p, slab_polygon(xs, th; position=pos, orientation=:vertical,
                                      side=side, justification=just, scale=sc, normalize=nrm);
                      color=(scol, salpha))
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

# (child recipe name => attribute defaults it presets over SlabInterval)
const _CHILD_DEFAULTS = Dict(
    :HalfEye        => (slab_type=:pdf,  side=:top,  show_slab=true,  show_interval=true,  show_point=true),
    :Eye            => (slab_type=:pdf,  side=:both, show_slab=true,  show_interval=true,  show_point=true),
    :CcdfInterval   => (slab_type=:ccdf, side=:top,  show_slab=true,  show_interval=true,  show_point=true),
    :CdfInterval    => (slab_type=:cdf,  side=:top,  show_slab=true,  show_interval=true,  show_point=true),
    :GradientInterval => (slab_type=:pdf, side=:both, show_slab=true, show_interval=true, show_point=true),
    :HistInterval   => (slab_type=:histogram, side=:top, show_slab=true, show_interval=true, show_point=true),
    :Slab           => (show_slab=true,  show_interval=false, show_point=false),
    :Interval       => (show_slab=false, show_interval=true,  show_point=false),
    :PointInterval  => (show_slab=false, show_interval=true,  show_point=true),
    :Spike          => (slab_type=:pdf,  show_slab=true, show_interval=false, show_point=false),
)

# Generate each child as a recipe that presets defaults then delegates to
# slabinterval!. This is ggdist's own "child = parent + defaults" mechanism.
for (T, defs) in _CHILD_DEFAULTS
    fname = Symbol(lowercase(String(T)))
    fname!  = Symbol(fname, :!)
    @eval begin
        @recipe($T, data) do scene
            merge(Attributes(; $(defs)...), Makie.default_theme(scene, SlabInterval))
        end
        function Makie.plot!(p::$T)
            slabinterval!(p, p.data;
                slab_type=p.slab_type, interval=p.interval, point=p.point, widths=p.widths,
                side=p.side, justification=p.justification, scale=p.scale, normalize=p.normalize,
                show_slab=p.show_slab, show_interval=p.show_interval, show_point=p.show_point,
                trim=p.trim, n=p.n, color=p.color, slab_color=p.slab_color, slab_alpha=p.slab_alpha,
                colormap=p.colormap, colorrange=p.colorrange,
                interval_linewidth=p.interval_linewidth, point_size=p.point_size)
            return p
        end
        # each public recipe type gets its own convert_arguments (Makie recipe
        # types do not share a supertype), all delegating to the shared converter
        Makie.convert_arguments(::Type{<:$T}, x::Distributions.UnivariateDistribution) = (_to_dist_args(x),)
        Makie.convert_arguments(::Type{<:$T}, x::AbstractVector{<:Real}) = (_to_dist_args(x),)
        Makie.convert_arguments(::Type{<:$T}, xs::AbstractVector{<:Distributions.UnivariateDistribution}) = (_to_dist_args(xs),)
    end
end

# Pre-summarised (ggdist geom_ path): draw supplied point/lower/upper directly.
function pointinterval(positions::AbstractVector{<:Real}, values::AbstractVector{<:Real},
                       lowers::AbstractVector{<:Real}, uppers::AbstractVector{<:Real}; kwargs...)
    fig = Figure()
    ax = Axis(fig[1,1])
    for i in eachindex(positions)
        a, b = interval_segment(lowers[i], uppers[i]; position=positions[i], orientation=:vertical)
        linesegments!(ax, [a, b]; linewidth=6, color=:black)
        scatter!(ax, [point_marker(values[i]; position=positions[i], orientation=:vertical)];
                 markersize=10, color=:black)
    end
    return Makie.FigureAxisPlot(fig, ax, ax.scene.plots[end])
end

# A slab needs a density; pre-summarised data has none — reject loudly.
_reject_slab_summary() = throw(ArgumentError(
    "slab needs a density; pre-summarised point/interval data has none. Pass samples or a Distribution, or use pointinterval for pre-summarised data."))
slab(::AbstractVector{<:Real}, ::AbstractVector{<:Real}, ::AbstractVector{<:Real}, ::AbstractVector{<:Real}) =
    _reject_slab_summary()
