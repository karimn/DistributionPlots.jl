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
