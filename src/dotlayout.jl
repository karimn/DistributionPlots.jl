"""
    dot_layout(d; ndots=50, binwidth=nothing)

Quantile dotplot layout (Wilkinson dot-density binning). Draws `ndots` evenly
spaced quantiles, then greedily packs them left-to-right into bins of width
`binwidth` (auto: the largest width that keeps every stack within range).
Returns `(x, y, binwidth)` where `y[i]` is the 1-based stack height of dot `i`.
"""
function dot_layout(d; ndots::Int=50, binwidth::Union{Nothing,Real}=nothing)
    dist = asdist(d)
    ps = (collect(1:ndots) .- 0.5) ./ ndots         # ndots representative quantiles
    vals = sort([quantile_at(dist, p) for p in ps])
    bw = binwidth === nothing ? _auto_binwidth(vals, ndots) : float(binwidth)

    x = Float64[]
    y = Int[]
    i = 1
    n = length(vals)
    while i ≤ n
        bin_start = vals[i]
        j = i
        while j ≤ n && vals[j] < bin_start + bw
            j += 1
        end
        members = vals[i:j-1]
        center = sum(members) / length(members)
        for (h, _) in enumerate(members)
            push!(x, center)
            push!(y, h)
        end
        i = j
    end
    return (x=x, y=y, binwidth=bw)
end

# Wilkinson's rule of thumb: binwidth from the data range and dot count.
function _auto_binwidth(vals::Vector{Float64}, ndots::Int)
    rng = vals[end] - vals[1]
    rng ≤ 0 && return 1.0
    return rng / (ndots / 2)      # ~2 dots per bin on average as a starting point
end
