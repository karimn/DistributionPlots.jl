using Makie

@recipe(Dots, data) do scene
    Attributes(
        ndots = 50,
        color = :black,
        scale = 0.9,
        markersize = 6,
        dodge = Makie.automatic,
        n_dodge = Makie.automatic,
        dodge_gap = 0.03,
        width = 1.0,
    )
end

# p.data already resolves to the (positions, dists) tuple produced by
# convert_arguments (Makie unwraps the recipe's single positional arg) —
# mirrors the pattern established in recipes/slabinterval.jl.
function Makie.plot!(p::Dots)
    args = p.data

    lift(args, p.ndots, p.scale, p.color, p.markersize,
         p.dodge, p.n_dodge, p.dodge_gap, p.width) do (positions, dists), nd, sc, col, msz,
            dg, ndg, dgap, wdt
        shift, shrink = dodge_placement(_or_nothing(dg), _or_nothing(ndg), dgap, wdt)
        sc = sc * shrink
        for (pos0, d) in zip(positions, dists)
            pos = pos0 + shift
            lay = dot_layout(d; ndots=nd)
            pts = [Point2f(pos + (y - 1) * lay.binwidth * sc, x) for (x, y) in zip(lay.x, lay.y)]
            scatter!(p, pts; markersize=msz, color=col)
        end
    end
    return p
end

Makie.convert_arguments(::Type{<:Dots}, x::Distributions.UnivariateDistribution) = (_to_dist_args(x),)
Makie.convert_arguments(::Type{<:Dots}, x::AbstractVector{<:Real}) = (_to_dist_args(x),)
Makie.convert_arguments(::Type{<:Dots}, xs::AbstractVector{<:Distributions.UnivariateDistribution}) =
    (_to_dist_args(xs),)
Makie.convert_arguments(::Type{<:Dots}, x::AbstractVector{<:Real}, y::AbstractVector{<:Real}) =
    (_to_dist_args(x, y),)

@recipe(DotsInterval, data) do scene
    Attributes(
        ndots = 50,
        show_interval = true,
        show_point = true,
        widths = [0.66, 0.95],
        point = :median,
        interval = :qi,
        color = :black,
        scale = 0.9,
        markersize = 6,
        interval_linewidth = Makie.automatic,
        point_size = 16,
        point_strokewidth = 1.5,
        point_strokecolor = :white,
        dodge = Makie.automatic,
        n_dodge = Makie.automatic,
        dodge_gap = 0.03,
        width = 1.0,
    )
end

function Makie.plot!(p::DotsInterval)
    args = p.data

    lift(args, p.ndots, p.widths, p.point, p.interval, p.show_interval, p.show_point,
         p.scale, p.color, p.markersize, p.interval_linewidth, p.point_size,
         p.point_strokewidth, p.point_strokecolor,
         p.dodge, p.n_dodge, p.dodge_gap, p.width) do (positions, dists),
            nd, ws, pt, iv, sint, spoint, sc, col, msz, ilw, psz, pstrokew, pstrokec,
            dg, ndg, dgap, wdt

        shift, shrink = dodge_placement(_or_nothing(dg), _or_nothing(ndg), dgap, wdt)
        sc = sc * shrink
        lws = interval_linewidths(ws, ilw)
        for (pos0, d) in zip(positions, dists)
            pos = pos0 + shift
            lay = dot_layout(d; ndots=nd)
            pts = [Point2f(pos + (y - 1) * lay.binwidth * sc, x) for (x, y) in zip(lay.x, lay.y)]
            scatter!(p, pts; markersize=msz, color=col)

            rows = point_interval(d; widths=ws, point=pt, interval=iv)
            if sint
                # widest (thinnest) first, so narrower (thicker) intervals draw on top
                for r in sort(rows; by=row -> row.width, rev=true)
                    a, b = interval_segment(r.lower, r.upper; position=pos, orientation=:vertical)
                    linesegments!(p, [a, b]; linewidth=lws[r.width], color=col)
                end
            end
            if spoint && !isempty(rows)
                scatter!(p, [point_marker(rows[1].value; position=pos, orientation=:vertical)];
                         markersize=psz, color=col, strokewidth=pstrokew, strokecolor=pstrokec)
            end
        end
    end
    return p
end

Makie.convert_arguments(::Type{<:DotsInterval}, x::Distributions.UnivariateDistribution) = (_to_dist_args(x),)
Makie.convert_arguments(::Type{<:DotsInterval}, x::AbstractVector{<:Real}) = (_to_dist_args(x),)
Makie.convert_arguments(::Type{<:DotsInterval}, xs::AbstractVector{<:Distributions.UnivariateDistribution}) =
    (_to_dist_args(xs),)
Makie.convert_arguments(::Type{<:DotsInterval}, x::AbstractVector{<:Real}, y::AbstractVector{<:Real}) =
    (_to_dist_args(x, y),)
