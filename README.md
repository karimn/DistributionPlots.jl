# DistributionPlots.jl

[![CI](https://github.com/karimn/DistributionPlots.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/karimn/DistributionPlots.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Makie recipes that mimic R's [`ggdist`](https://mjskay.github.io/ggdist/) —
slab, interval, pointinterval, lineribbon, and the dots family — for visualizing
distributions and uncertainty. Plots [`RVars.jl`](https://github.com/karimn/RVars.jl)
`RVar` samples, `Distributions.jl` distributions, `MCMCChains.Chains`, and
raw sample vectors through one interface. **Built on Makie** (not Plots.jl).

## Install

```julia
using Pkg
Pkg.add(url="https://github.com/karimn/RVars.jl")
Pkg.add(url="https://github.com/karimn/DistributionPlots.jl")
```

## Quick start

```julia
using DistributionPlots, CairoMakie, RVars

rv = RVar(randn(1000, 3))     # 3 parameters, 1000 draws each
halfeye(rv)                          # density slab + interval + point, per parameter
pointinterval(rv)                    # just the point + nested intervals
```

You can also pass raw draws plus the position to group them by — ggdist's
`stat_halfeye(aes(x = arm, y = value))`, one distribution per distinct `x`:

```julia
halfeye(x, values)
```

Sample vectors and analytic distributions work too:

```julia
using Distributions
halfeye(randn(2000))
halfeye(Normal(0, 1))
```

`point_interval` is a standalone tidybayes-style summary (no plotting):

```julia
point_interval(randn(2000); widths=[0.66, 0.95], point=:median, interval=:hdi)
# Vector of NamedTuples: (value, lower, upper, width, point, interval)
```

## Named dimensions and categories

A chain records only `a[1,2]`; that axis 1 is trials and axis 2 is arms lives in
the model. [`RVars.jl`](https://github.com/karimn/RVars.jl) lets you declare that
when parameters are extracted, and DistributionPlots uses it for tick labels,
colours, grouping and facets.

```julia
p = RVar(chn; dims   = (a = (:trial, :arm),),
              labels = (arm = ["control", "drug"],))
```

### One axis

A single Makie axis can show one dimension, so pass `dim` to say which:

```julia
fig, ax, plt = halfeye(p.a; dim = :arm)   # the other dimensions are pooled
ax.xticks = dimticks(p.a, :arm)           # ("control", "drug")
```

`dimticks(x, dim)` takes labels from `RVars.dimlabels`, falling back to
`variables(x)` for a named vector random variable and to integer positions
otherwise. Ticks must be applied by the caller because a Makie recipe draws into
a `Scene` and cannot reach the enclosing `Axis`.

A rank-2 or higher `RVar` without `dim` is an error naming the dimensions it
found, rather than a silently mis-flattened plot.

### Colours, grouping and facets

Faceting and legends belong to the grammar layer, not to a recipe — so go
through `RVars.gather_draws`, which flattens a labelled `RVar` into a long table
whose columns are named after the dimensions and hold the categories themselves:

```julia
using AlgebraOfGraphics

tbl = RVars.gather_draws(p.a)   # (trial, arm, chain, draw, value)

data(tbl) * mapping(:arm, :value; color = :arm, layout = :trial) *
  visual(HalfEye) |> draw
```

Every facet, colour, legend entry, axis label and tick label there comes from the
dimension names and labels. `mapping` accepts the usual AoG aesthetics —
`color`, `layout`, `row`, `col`, `dodge`.

ggdist's colour/fill split is spelled `color`/`slab_color` here: `color` drives
the point and interval, `slab_color` the density, and the two map independently.
So ggdist's `aes(fill = arm)` is `mapping(slab_color = :arm)`. The name follows
the recipe's own `slab_type`/`slab_alpha`/`show_slab` vocabulary and Makie's
`color`/`*color` convention — neither Makie nor AlgebraOfGraphics has a `fill`
aesthetic.

Dodging needs `n_dodge` passed explicitly, because AoG calls a recipe once per
group with a scalar `dodge` and never reports how many groups there are:

```julia
data(tbl) * mapping(:trial, :value; color = :arm, dodge = :arm) *
  visual(HalfEye; n_dodge = 2) |> draw
```

## AlgebraOfGraphics

Loading `AlgebraOfGraphics` alongside `DistributionPlots` activates a package
extension registering, for every public recipe type, the `aesthetic_mapping` AoG
uses to place arguments and attributes on aesthetics, and the `legend_elements`
it needs to build a legend. Recipes take `(x, values)`, matching AoG's usual
two-column `mapping(:x, :y)` form; the single-column `mapping(:y) *
visual(Interval)` form also works and draws one distribution.

Tested against AlgebraOfGraphics v0.13.1; its `aesthetic_mapping` API is not yet
fully stable across releases, so pin or check compatibility if you hit
`MethodError`/`aesthetic_mapping` errors on a different version.

## Development

`RVars.jl` is unregistered, so the test environment develops it (and this package)
explicitly rather than using `Pkg.test()`. To run the suite locally:

```bash
julia --project=test -e 'using Pkg; Pkg.add(url="https://github.com/karimn/RVars.jl"); Pkg.develop(path="."); Pkg.instantiate()'
julia --project=test test/runtests.jl
```

(CI runs these same two steps — see `.github/workflows/CI.yml`.)

## Notes

- Chains support requires `using RVars` alongside `using MCMCChains`.
- `halfeye(chn)` plots one distribution per scalar parameter. To get shaped
  parameters with named dimensions, extract them first with
  `rvars(chn; dims = ...)` and plot the resulting `RVar`.
- Requires RVars 0.5.1+ (dimension names and labels, `gather_draws`).
- Requires Julia 1.10+.
