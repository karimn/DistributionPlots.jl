using Makie

@recipe(Dots, data) do scene
    Attributes(
        ndots = 50,
        color = :black,
        scale = 0.9,
        markersize = 6,
    )
end

# p.data already resolves to the (positions, dists) tuple produced by
# convert_arguments (Makie unwraps the recipe's single positional arg) —
# mirrors the pattern established in recipes/slabinterval.jl.
function Makie.plot!(p::Dots)
    args = p.data

    lift(args, p.ndots, p.scale, p.color, p.markersize) do (positions, dists), nd, sc, col, msz
        for (pos, d) in zip(positions, dists)
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
        interval_linewidth = 6,
        point_size = 10,
    )
end

function Makie.plot!(p::DotsInterval)
    args = p.data

    lift(args, p.ndots, p.widths, p.point, p.interval, p.show_interval, p.show_point,
         p.scale, p.color, p.markersize, p.interval_linewidth, p.point_size) do (positions, dists),
            nd, ws, pt, iv, sint, spoint, sc, col, msz, ilw, psz

        for (pos, d) in zip(positions, dists)
            lay = dot_layout(d; ndots=nd)
            pts = [Point2f(pos + (y - 1) * lay.binwidth * sc, x) for (x, y) in zip(lay.x, lay.y)]
            scatter!(p, pts; markersize=msz, color=col)

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

Makie.convert_arguments(::Type{<:DotsInterval}, x::Distributions.UnivariateDistribution) = (_to_dist_args(x),)
Makie.convert_arguments(::Type{<:DotsInterval}, x::AbstractVector{<:Real}) = (_to_dist_args(x),)
Makie.convert_arguments(::Type{<:DotsInterval}, xs::AbstractVector{<:Distributions.UnivariateDistribution}) =
    (_to_dist_args(xs),)
