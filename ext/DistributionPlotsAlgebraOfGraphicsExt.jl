module DistributionPlotsAlgebraOfGraphicsExt

using DistributionPlots
using AlgebraOfGraphics
using Makie

# AlgebraOfGraphics can already wrap any Makie recipe via visual(RecipeType), but it
# needs to know which aesthetic (x/y/...) each positional argument maps to. Our
# recipes take (positions, distributions) as (x, y)-like continuous data, so register
# the straightforward 1 => AesX, 2 => AesY mapping for each public recipe type.
# Without this, AoG raises "No aesthetic mapping defined yet for plot type ...".
# (AlgebraOfGraphics v0.13.1's `aesthetic_mapping` splats one `ScientificType` per
# positional plot argument, so the fallback must accept a variable number of them.)
import AlgebraOfGraphics: aesthetic_mapping
using AlgebraOfGraphics: AesX, AesY, ScientificType, dictionary
for T in (DistributionPlots.SlabInterval, DistributionPlots.HalfEye,
          DistributionPlots.Interval, DistributionPlots.PointInterval)
    @eval aesthetic_mapping(::Type{<:$T}, ::ScientificType...) = dictionary([1 => AesX, 2 => AesY])
end

end # module
