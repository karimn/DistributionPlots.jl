using Makie

# Shared converter: normalise any input to (positions, Vector{AbstractDist}).
# A single distribution/sample-vector is placed at position 1. Extensions
# (RVars, MCMCChains) add methods that expand multi-element inputs to 1:k.
_to_dist_args(x::Distributions.UnivariateDistribution) = ([1.0], AbstractDist[asdist(x)])
_to_dist_args(x::AbstractVector{<:Real}) = ([1.0], AbstractDist[asdist(x)])
function _to_dist_args(xs::AbstractVector{<:Distributions.UnivariateDistribution})
    return (collect(1.0:length(xs)), AbstractDist[asdist(x) for x in xs])
end

# Two-argument (x, values) form: raw draws grouped by their position, one
# distribution per distinct x. This is ggdist's `stat_halfeye(aes(x=arm, y=value))`
# and the shape AlgebraOfGraphics feeds a recipe — it maps a categorical column to
# integer codes, sets the tick labels on the Axis itself, and hands the recipe every
# row of a facet at once, leaving the grouping to us.
function _to_dist_args(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    length(x) == length(y) || throw(DimensionMismatch(
        "x and values must have equal length, got $(length(x)) and $(length(y))"))
    positions = sort!(unique(x))
    slot = Dict(p => i for (i, p) in enumerate(positions))
    groups = [Float64[] for _ in positions]
    for (xi, yi) in zip(x, y)
        push!(groups[slot[xi]], yi)
    end
    return (collect(Float64, positions), AbstractDist[asdist(g) for g in groups])
end

# Makie signals "not set" with `automatic`; the geometry helpers take `nothing`.
_or_nothing(x) = x === Makie.automatic ? nothing : x

@recipe(SlabInterval, data) do scene
    Attributes(
        slab_type = :pdf,          # :pdf, :cdf, :ccdf, :histogram
        interval = :qi,            # :qi, :hdci, :hdi
        point = :median,           # :mean, :median, :mode
        widths = [0.66, 0.95],
        side = :top,               # :top, :bottom, :both
        justification = 0.0,
        scale = 0.9,
        normalize = :all,          # :all (one global max across slabs), :each (per-slab), :none
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
        # Side-by-side groups within one position slot. AlgebraOfGraphics drives
        # these via the AesDodgeX aesthetic; see dodge_placement in geometry.jl for
        # why n_dodge has to be supplied rather than inferred.
        dodge = Makie.automatic,
        n_dodge = Makie.automatic,
        dodge_gap = 0.03,
        width = 1.0,
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
Makie.convert_arguments(::Type{<:SlabInterval}, x::AbstractVector{<:Real}, y::AbstractVector{<:Real}) =
    (_to_dist_args(x, y),)

function Makie.plot!(p::SlabInterval)
    # p.data already resolves to the (positions, dists) tuple produced by
    # convert_arguments (Makie unwraps the recipe's single positional arg).
    args = p.data

    lift(args, p.slab_type, p.interval, p.point, p.widths, p.side, p.justification,
         p.scale, p.normalize, p.trim, p.n, p.show_slab, p.show_interval, p.show_point,
         p.color, p.slab_color, p.slab_alpha, p.interval_linewidth, p.point_size,
         p.dodge, p.n_dodge, p.dodge_gap, p.width) do (positions, dists), st, iv, pt,
            ws, side, just, sc, nrm, tr, n, sslab, sint, spoint, col, scol, salpha, ilw, psz,
            dg, ndg, dgap, wdt

        # One placement for the whole plot: AoG calls the recipe once per dodge group.
        # The shift applies to the position rather than to `justification`, which
        # offsets only the slab — a dodged group moves slab, interval and point
        # together — and the slab narrows to its share of the slot.
        shift, shrink = dodge_placement(_or_nothing(dg), _or_nothing(ndg), dgap, wdt)
        sc = sc * shrink

        # Precompute slab curves so `normalize=:all` can divide every slab by ONE
        # global maximum thickness (ggdist's "all"), preserving relative peak heights
        # across slabs. `:each` normalizes each slab to its own max; `:none` leaves raw.
        curves = sslab ? [slab_curve(d; kind=st, n=n, trim=tr) for d in dists] : nothing
        gmax = 0.0
        if sslab
            for (_, th) in curves
                isempty(th) || (gmax = max(gmax, maximum(th)))
            end
        end
        for (i, (pos0, d)) in enumerate(zip(positions, dists))
            pos = pos0 + shift
            # slab drawn first so the point/interval below render on top of it
            if sslab
                xs, th = curves[i]
                poly!(p, slab_polygon(xs, th; position=pos, orientation=:vertical,
                                      side=side, justification=just, scale=sc,
                                      normalize=nrm, globalmax=gmax);
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
                interval_linewidth=p.interval_linewidth, point_size=p.point_size,
                dodge=p.dodge, n_dodge=p.n_dodge, dodge_gap=p.dodge_gap, width=p.width)
            return p
        end
        # each public recipe type gets its own convert_arguments (Makie recipe
        # types do not share a supertype), all delegating to the shared converter
        Makie.convert_arguments(::Type{<:$T}, x::Distributions.UnivariateDistribution) = (_to_dist_args(x),)
        Makie.convert_arguments(::Type{<:$T}, x::AbstractVector{<:Real}) = (_to_dist_args(x),)
        Makie.convert_arguments(::Type{<:$T}, xs::AbstractVector{<:Distributions.UnivariateDistribution}) = (_to_dist_args(xs),)
        Makie.convert_arguments(::Type{<:$T}, x::AbstractVector{<:Real}, y::AbstractVector{<:Real}) = (_to_dist_args(x, y),)
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
