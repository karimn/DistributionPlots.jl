module DistributionPlotsMCMCChainsExt

using DistributionPlots
using DistributionPlots: _to_dist_args
using MCMCChains: Chains
using RVars: RVar
using Makie

# Route Chains through RVars' own Chains conversion (single tested path), then reuse
# the RVars _to_dist_args method.
#
# RVars 0.4 made `RVar(chn)` return a NamedTuple of shaped, per-parameter RVars rather
# than one flat named vector. Plotting a whole fit means "one distribution per scalar
# parameter", which is the flat view, so ask for it explicitly. Extracting shaped
# parameters (and their dimension names and labels) is a separate, deliberate step the
# user makes with rvars(chn; dims=...) — plot the resulting RVar directly.
DistributionPlots._to_dist_args(chn::Chains) = _to_dist_args(RVar(chn; flat = true))

const _RECIPES = (SlabInterval, HalfEye, Eye, CcdfInterval, CdfInterval,
                  GradientInterval, HistInterval, Slab, Interval, PointInterval,
                  Spike, Dots, DotsInterval)

for T in _RECIPES
    @eval Makie.convert_arguments(::Type{<:$T}, chn::Chains) = (_to_dist_args(chn),)
end

end # module
