using Makie

@recipe(LineRibbon, x, dists) do scene
    Attributes(
        widths = [0.5, 0.8, 0.95],
        point = :median,
        interval = :qi,
        colormap = :blues,
        color = :black,
        linewidth = 2,
    )
end

function Makie.plot!(p::LineRibbon)
    lift(p.x, p.dists, p.widths, p.point, p.interval, p.colormap, p.color,
         p.linewidth) do xs, dists, ws, pt, iv, cmap, col, lw
        wsorted = sort(ws; rev=true)              # widest band drawn first (back), narrower on top
        ncol = length(wsorted)
        colors = Makie.cgrad(cmap, max(ncol, 2))
        for (k, w) in enumerate(wsorted)
            los = Float64[]
            his = Float64[]
            for d in dists
                r = point_interval(d; widths=[w], point=pt, interval=iv)[1]
                push!(los, r.lower)
                push!(his, r.upper)
            end
            # narrower bands (higher k) are drawn later (on top) and get a
            # progressively more saturated colour so they read as "more
            # certain" nested inside the wider, fainter bands behind them
            band!(p, xs, los, his; color=colors[k])
        end
        centers = [point_interval(d; widths=[first(wsorted)], point=pt, interval=iv)[1].value for d in dists]
        lines!(p, xs, centers; color=col, linewidth=lw)
    end
    return p
end

# predictor + per-x sample vectors → wrap each entry as an AbstractDist
Makie.convert_arguments(::Type{<:LineRibbon}, x::AbstractVector{<:Real},
                        dists::AbstractVector{<:AbstractVector{<:Real}}) =
    (collect(Float64, x), AbstractDist[asdist(d) for d in dists])
Makie.convert_arguments(::Type{<:LineRibbon}, x::AbstractVector{<:Real},
                        dists::AbstractVector{<:Distributions.UnivariateDistribution}) =
    (collect(Float64, x), AbstractDist[asdist(d) for d in dists])
