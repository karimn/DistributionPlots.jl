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

function interval_segment(lower::Real, upper::Real; position::Real, orientation::Symbol=:vertical)
    o = Val(orientation)
    return (_pt(lower, position, o), _pt(upper, position, o))
end

point_marker(value::Real; position::Real, orientation::Symbol=:vertical) =
    _pt(value, position, Val(orientation))
