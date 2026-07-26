using Makie: Point2f

function normalize_thickness(thickness::AbstractVector{<:Real}, mode::Symbol;
                             globalmax::Union{Nothing,Real}=nothing)
    if mode === :none
        return collect(Float64, thickness)
    elseif mode === :each
        m = isempty(thickness) ? 0.0 : maximum(thickness)
        return m > 0 ? collect(Float64, thickness ./ m) : collect(Float64, thickness)
    elseif mode === :all
        # Divide by ONE global maximum when the recipe supplies it (across ALL slabs,
        # preserving relative peak heights — ggdist's "all"). When called standalone
        # with no global max, fall back to this vector's own max.
        m = globalmax === nothing ? (isempty(thickness) ? 0.0 : maximum(thickness)) : float(globalmax)
        return m > 0 ? collect(Float64, thickness ./ m) : collect(Float64, thickness)
    else
        throw(ArgumentError("normalize_thickness: unknown mode :$mode (expected :all, :each, or :none)"))
    end
end

# Map a (value, offset) pair to a Point2f under an orientation.
# `along` is the value-axis coordinate; `perp` is the thickness/position axis.
_pt(along, perp, ::Val{:vertical}) = Point2f(perp, along)
_pt(along, perp, ::Val{:horizontal}) = Point2f(along, perp)
_pt(along, perp, ::Val{o}) where {o} = throw(ArgumentError(
    "orientation must be :vertical or :horizontal, got :$o"))

function slab_polygon(xs::AbstractVector{<:Real}, thickness::AbstractVector{<:Real};
                      position::Real, orientation::Symbol=:vertical,
                      side::Symbol=:top, justification::Real=0.0,
                      scale::Real=1.0, normalize::Symbol=:all, globalmax=nothing)
    th = normalize_thickness(thickness, normalize; globalmax=globalmax) .* scale
    o = Val(orientation)
    base = position + justification
    top = [_pt(xs[i], base + _sideoffset(th[i], side, :hi), o) for i in eachindex(xs)]
    bot = [_pt(xs[i], base + _sideoffset(th[i], side, :lo), o) for i in reverse(eachindex(xs))]
    return vcat(top, bot)
end

_sideoffset(t, ::Val, ::Any) = t
function _sideoffset(t, side::Symbol, which::Symbol)
    if side === :top
        return which === :hi ? t : 0.0
    elseif side === :bottom
        return which === :hi ? 0.0 : -t
    elseif side === :both
        return which === :hi ? t/2 : -t/2
    else
        throw(ArgumentError("slab: unknown side :$side (expected :top, :bottom, or :both)"))
    end
end

# Side-by-side placement within one position slot, matching Makie's barplot so a
# dodged slab lines up with a dodged bar drawn on the same axis (same two formulas,
# inlined rather than called because they are Makie internals).
#
# Returns (shift, shrink): how far to move this group's centre, and the factor its
# width must shrink by so the groups sit side by side instead of on top of one
# another — dodged bars narrow the same way.
dodge_placement(::Nothing, n_dodge, dodge_gap, width) = (0.0, 1.0)
function dodge_placement(dodge::Integer, n_dodge::Integer, dodge_gap, width)
    n_dodge <= 1 && return (0.0, 1.0)
    1 <= dodge <= n_dodge || throw(ArgumentError(
        "dodge=$dodge is out of range for n_dodge=$n_dodge"))
    dodge_width = (1 - (n_dodge - 1) * dodge_gap) / n_dodge
    shift = ((dodge_width - 1) / 2 + (dodge - 1) * (dodge_width + dodge_gap)) * width
    return (shift, dodge_width)
end
# `n_dodge` cannot be inferred when AlgebraOfGraphics drives the plot: it splits the
# data by group and calls the recipe once per group, passing a *scalar* `dodge` and
# leaving `n_dodge` automatic, so no single call can see how many groups there are.
# Guessing would silently mis-size every group, so ask for it instead.
dodge_placement(dodge::Integer, ::Nothing, dodge_gap, width) = throw(ArgumentError(
    "dodge=$dodge was given without n_dodge, and the number of dodge groups cannot " *
    "be inferred from a single group index. Pass n_dodge explicitly (e.g. " *
    "visual(HalfEye; n_dodge=2) under AlgebraOfGraphics)."))
# One plot call draws one dodge group, so a per-element dodge vector has no meaning
# here — unlike barplot, where a single call draws many independently dodged bars.
dodge_placement(dodge, n_dodge, dodge_gap, width) = throw(ArgumentError(
    "dodge must be a single group index, got $(typeof(dodge)). Each plot call draws " *
    "one dodge group; call the recipe once per group (which is what " *
    "AlgebraOfGraphics does)."))

function interval_segment(lower::Real, upper::Real; position::Real, orientation::Symbol=:vertical)
    o = Val(orientation)
    return (_pt(lower, position, o), _pt(upper, position, o))
end

point_marker(value::Real; position::Real, orientation::Symbol=:vertical) =
    _pt(value, position, Val(orientation))

# Per-width line weight for nested credible intervals, plus the draw order they
# belong in. ggdist's nested look comes from two things at once: the widest
# (least certain) interval is thin and drawn first, so narrower ones layer on
# top of it instead of being hidden inside it.
#
# `interval_linewidth` can be:
#   - `Makie.automatic`: derive a thick-to-thin ramp, thickest for the
#     narrowest width, thinnest for the widest.
#   - a vector: one entry per distinct width, given in the same
#     widest-to-narrowest order the intervals are drawn in.
#   - a scalar: every width gets that same line weight (pre-#5 behaviour).
function interval_linewidths(widths::AbstractVector{<:Real}, interval_linewidth)
    sorted = sort(unique(Float64.(widths)); rev=true)   # widest (drawn first/back) → narrowest (last/front)
    k = length(sorted)
    if interval_linewidth isa AbstractVector
        length(interval_linewidth) == k || throw(ArgumentError(
            "interval_linewidth vector must have one entry per distinct width " *
            "($k here), got $(length(interval_linewidth))"))
        lws = Float64.(interval_linewidth)
    elseif interval_linewidth === Makie.automatic
        thick, thin = 8.0, 3.0
        lws = k == 1 ? [(thick + thin) / 2] :
              [thin + (thick - thin) * (i - 1) / (k - 1) for i in 1:k]
    else
        lws = fill(Float64(interval_linewidth), k)
    end
    return Dict(sorted[i] => lws[i] for i in 1:k)
end
