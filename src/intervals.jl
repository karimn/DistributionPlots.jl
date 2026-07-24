const PIRow = @NamedTuple{value::Float64, lower::Float64, upper::Float64,
                          width::Float64, point::Symbol, interval::Symbol}

function point_summary(d::AbstractDist, point::Symbol)
    if point === :mean
        return _mean(d)
    elseif point === :median
        return quantile_at(d, 0.5)
    elseif point === :mode
        return _mode(d)           # defined in Task 5
    else
        throw(ArgumentError("point_interval: unknown point summary :$point (expected :mean, :median, or :mode)"))
    end
end

_mean(d::SampleDist) = mean(d.samples)
_mean(d::AnalyticDist) = mean(d.dist)

"""
    point_interval(d; widths=[0.66, 0.95], point=:median, interval=:qi)

Compute a point summary and one or more probability intervals, returned as a
`Vector{NamedTuple}` (Tables.jl row table) with columns
`value, lower, upper, width, point, interval`.
"""
function point_interval(d; widths=[0.66, 0.95], point::Symbol=:median, interval::Symbol=:qi)
    dist = asdist(d)
    val = point_summary(dist, point)
    rows = PIRow[]
    for w in widths
        (0 < w < 1) || throw(ArgumentError("point_interval: width must be in (0,1), got $w"))
        for (lo, hi) in _interval_pieces(dist, w, interval)
            push!(rows, (value=val, lower=lo, upper=hi, width=float(w),
                         point=point, interval=interval))
        end
    end
    return rows
end

# Returns a vector of (lower, upper) pieces. qi/hdci → one piece.
function _interval_pieces(d::AbstractDist, w::Real, interval::Symbol)
    if interval === :qi
        α = (1 - w) / 2
        return [(quantile_at(d, α), quantile_at(d, 1 - α))]
    elseif interval === :hdci
        return [_hdci(d, w)]
    elseif interval === :hdi
        return _hdi(d, w)          # defined in Task 5
    else
        throw(ArgumentError("point_interval: unknown interval :$interval (expected :qi, :hdci, or :hdi)"))
    end
end

# Highest-density *continuous* interval: narrowest contiguous window holding mass w.
function _hdci(d::SampleDist, w::Real)
    x = sort(d.samples)
    n = length(x)
    k = max(1, round(Int, w * n))
    k ≥ n && return (x[1], x[end])
    best_lo, best_hi, best_span = x[1], x[k], x[k] - x[1]
    for i in 1:(n - k)
        span = x[i + k] - x[i]
        if span < best_span
            best_span, best_lo, best_hi = span, x[i], x[i + k]
        end
    end
    return (best_lo, best_hi)
end

# Analytic hdci: minimise (quantile(u+w) - quantile(u)) over u in [0, 1-w].
function _hdci(d::AnalyticDist, w::Real)
    us = range(0.0, 1 - w; length=1001)
    best_u, best_span = 0.0, Inf
    for u in us
        span = quantile_at(d, u + w) - quantile_at(d, u)
        if span < best_span
            best_span, best_u = span, u
        end
    end
    return (quantile_at(d, best_u), quantile_at(d, best_u + w))
end
