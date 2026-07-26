module DistributionPlotsRVarsExt

using DistributionPlots
using DistributionPlots: _to_dist_args, AbstractDist, asdist
using RVars
using RVars: RVar, draws, variables, dimnames, dimlabels
using Makie
using Distributions

# Resolve a dimension given by name (or position) to an axis index, mirroring the
# errors RVars' own (unexported) _dim_index raises.
function _dim_index(x::RVar, dim::Symbol)
    dn = dimnames(x)
    dn === nothing && throw(ArgumentError(
        "dim=:$dim was given but this RVar has no dimension names. Attach them at " *
        "extraction with rvars(chn; dims=(a=(:trial, :arm),)), or index by position."))
    d = findfirst(isequal(dim), dn)
    d === nothing && throw(ArgumentError(
        "unknown dimension :$dim; available: $(join(string.(dn), ", "))"))
    return d
end
function _dim_index(x::RVar, dim::Integer)
    1 <= dim <= ndims(x) || throw(ArgumentError(
        "dim=$dim is out of range for a rank-$(ndims(x)) RVar"))
    return Int(dim)
end

"""
    dimticks(x::RVar, dim = 1)

Makie tick specification `(positions, labels)` for the axis `x` is plotted against,
taking the labels from [`RVars.dimlabels`](@ref) when the dimension is labelled, from
`variables(x)` for a named vector random variable, and from the integer positions
otherwise. `dim` may be a dimension name or a position.

A Makie recipe draws into a `Scene` and cannot reach the enclosing `Axis`, so tick
labels have to be applied by the caller:

```julia
fig, ax, plt = halfeye(rv; dim = :arm)
ax.xticks = dimticks(rv, :arm)
```

`AlgebraOfGraphics` sets equivalent ticks by itself, so this is only needed on the
direct plotting path.
"""
function DistributionPlots.dimticks(x::RVar, dim = 1)
    ndims(x) == 0 && throw(ArgumentError(
        "a scalar (rank-0) RVar is drawn at a single position and has no dimension to label"))
    d = _dim_index(x, dim)
    k = size(x, d)
    labs = dimlabels(x, d)
    if labs === nothing && ndims(x) == 1
        labs = variables(x)   # element names on the flat/named vector path
    end
    return (1:k, labs === nothing ? string.(1:k) : string.(labs))
end

# `dim` selects which axis runs along the position axis; the remaining axes are
# pooled into each position's distribution. That is deliberately the same semantics
# as `data(gather_draws(x)) * mapping(:arm, :value)` under AlgebraOfGraphics, which
# also groups by x and pools everything else — so both paths draw the same picture.
function DistributionPlots._to_dist_args(x::RVar; dim = nothing)
    N = ndims(x)
    if N == 0
        dim === nothing || throw(ArgumentError(
            "dim=$(repr(dim)) was given for a scalar (rank-0) RVar, which has no dimensions"))
        return ([1.0], AbstractDist[asdist(vec(draws(x)))])
    end
    if N > 1 && dim === nothing
        dn = dimnames(x)
        named = dn === nothing ? "" : " (" * join(string.(dn), ", ") * ")"
        throw(ArgumentError(
            "this RVar has $N dimensions$named; a single axis can show only one. " *
            "Pass dim=:name to choose which runs along the axis (the others are " *
            "pooled), slice first (e.g. x[trial=1]), or use RVars.gather_draws with " *
            "AlgebraOfGraphics to map the other dimensions to colour or facets."))
    end
    d = dim === nothing ? 1 : _dim_index(x, dim)
    k = size(x, d)
    dists = AbstractDist[
        asdist(vec(draws(x[ntuple(j -> j == d ? i : Colon(), N)...])))
        for i in 1:k
    ]
    return (collect(1.0:k), dists)
end

# Register convert_arguments for every public recipe type. `used_attributes` is what
# routes `dim` into the converter: Makie strips those attributes from the plot and
# passes them to convert_arguments as keywords, which is the only way a converter can
# see an attribute at all.
const _RECIPES = (SlabInterval, HalfEye, Eye, CcdfInterval, CdfInterval,
                  GradientInterval, HistInterval, Slab, Interval, PointInterval,
                  Spike, Dots, DotsInterval)

for T in _RECIPES
    @eval begin
        Makie.used_attributes(::Type{<:$T}, ::RVar) = (:dim,)
        Makie.convert_arguments(::Type{<:$T}, x::RVar; dim = nothing) =
            (_to_dist_args(x; dim = dim),)
    end
end

end # module
