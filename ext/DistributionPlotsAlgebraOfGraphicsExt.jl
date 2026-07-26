module DistributionPlotsAlgebraOfGraphicsExt

using DistributionPlots
using AlgebraOfGraphics
using Makie

# AlgebraOfGraphics can wrap any Makie recipe via visual(RecipeType), but it needs two
# things from us: which aesthetic each argument and attribute maps to, and what a
# legend entry for the recipe looks like.
#
# Our recipes take (x, values) — raw draws plus the position to group them by, i.e.
# ggdist's stat_halfeye(aes(x = arm, y = value)). AoG maps a categorical x column to
# integer codes and sets the Axis tick labels itself, splits the table by the grouping
# aesthetics, and calls the recipe once per group with an already-resolved colour. So
# colours, grouping, faceting and the legend all come from the mapping — the recipe
# only ever draws one group into one axis.
#
# (AlgebraOfGraphics v0.13.1's `aesthetic_mapping` splats one `ScientificType` per
# positional plot argument, so the fallback must accept a variable number of them.)
#
# Issues #4/#8: our recipes swap which argument is drawn on which axis when
# `orientation = :horizontal`, but a static `1 => AesX, 2 => AesY` mapping doesn't
# follow — every label, tick and dodge direction stays as if it were still vertical.
# AoG's fix for exactly this (see Makie's own Violin/BoxPlot) is an attribute-dependent
# mapping: instead of an `Aes` value, map to `attrkey => dictionary(value => Aes...)`,
# and AoG resolves it against the plot's actual attribute value at draw time. That
# resolution looks the attribute up in the attributes AoG has for this plot — which
# only holds what was explicitly passed to `visual(...; kwargs...)`, not the recipe's
# own attribute defaults — so `mandatory_attributes` is needed to supply :vertical
# when the caller didn't set `orientation` at all.
import AlgebraOfGraphics: aesthetic_mapping, legend_elements, mandatory_attributes
using AlgebraOfGraphics: AesX, AesY, AesColor, AesDodgeX, AesDodgeY, ScientificType,
                         dictionary, MixedArguments, PolyElement, LineElement, MarkerElement

const _ORIENT_XY = :orientation => dictionary([:vertical => AesX, :horizontal => AesY])
const _ORIENT_YX = :orientation => dictionary([:vertical => AesY, :horizontal => AesX])
const _ORIENT_DODGE = :orientation => dictionary([:vertical => AesDodgeX, :horizontal => AesDodgeY])

# The slabinterval family inherits SlabInterval's attributes, so all of them have a
# `slab_color` distinct from `color`. Dots and DotsInterval draw no slab and have only
# `color` — declaring slab_color for them would offer a mapping the plot cannot accept.
const _SLAB_FAMILY = (DistributionPlots.SlabInterval, DistributionPlots.HalfEye,
                      DistributionPlots.Eye, DistributionPlots.CcdfInterval,
                      DistributionPlots.CdfInterval, DistributionPlots.GradientInterval,
                      DistributionPlots.HistInterval, DistributionPlots.Slab,
                      DistributionPlots.Spike, DistributionPlots.Interval,
                      DistributionPlots.PointInterval)
const _DOTS_FAMILY = (DistributionPlots.Dots, DistributionPlots.DotsInterval)

for T in _SLAB_FAMILY
    @eval aesthetic_mapping(::Type{<:$T}, ::ScientificType...) = dictionary([
        1 => _ORIENT_XY,
        2 => _ORIENT_YX,
        # `color` drives the point and interval, `slab_color` the density itself —
        # ggdist's colour/fill split, so mapping(color = :arm) tints the intervals and
        # mapping(slab_color = :arm) tints the slabs, independently.
        :color => AesColor,
        :slab_color => AesColor,
        # Native dodge: AoG hands the recipe a scalar `dodge` group index and leaves
        # `n_dodge` unset, so pass n_dodge yourself — visual(HalfEye; n_dodge = 2).
        :dodge => _ORIENT_DODGE,
    ])
    @eval mandatory_attributes(::Type{<:$T}) = dictionary([:orientation => :vertical])
end

for T in _DOTS_FAMILY
    @eval aesthetic_mapping(::Type{<:$T}, ::ScientificType...) = dictionary([
        1 => _ORIENT_XY, 2 => _ORIENT_YX, :color => AesColor, :dodge => _ORIENT_DODGE,
    ])
    @eval mandatory_attributes(::Type{<:$T}) = dictionary([:orientation => :vertical])
end

_legend_color(scale_args) = get(scale_args, :color, :black)

# A filled patch for anything that draws a slab; a line plus a marker for the
# interval-only recipes, which have no area to show; a marker for the dot families.
for T in (DistributionPlots.SlabInterval, DistributionPlots.HalfEye,
          DistributionPlots.Eye, DistributionPlots.CcdfInterval,
          DistributionPlots.CdfInterval, DistributionPlots.GradientInterval,
          DistributionPlots.HistInterval, DistributionPlots.Slab,
          DistributionPlots.Spike)
    @eval legend_elements(::Type{<:$T}, attributes, scale_args::MixedArguments) =
        [PolyElement(color = get(scale_args, :slab_color, _legend_color(scale_args)))]
end

for T in (DistributionPlots.Interval, DistributionPlots.PointInterval)
    @eval legend_elements(::Type{<:$T}, attributes, scale_args::MixedArguments) =
        [LineElement(color = _legend_color(scale_args)),
         MarkerElement(color = _legend_color(scale_args), marker = :circle)]
end

for T in _DOTS_FAMILY
    @eval legend_elements(::Type{<:$T}, attributes, scale_args::MixedArguments) =
        [MarkerElement(color = _legend_color(scale_args), marker = :circle)]
end

end # module
