# DistributionPlots.jl

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

## AlgebraOfGraphics

Loading `AlgebraOfGraphics` alongside `DistributionPlots` activates a package
extension so `AlgebraOfGraphics.visual` can wrap our recipe types (e.g.
`visual(Interval)`, `visual(HalfEye)`) inside an AoG spec. This is a minimal
integration point — it registers the `aesthetic_mapping` AoG needs to place a
recipe's positional arguments on the x/y axes — not a full set of AoG-native
analyses; those are out of scope for v1. Our recipes' `convert_arguments`
contract takes one positional argument (a distribution, a sample vector, or a
vector of distributions), so map a single data column (e.g.
`mapping(:y) * visual(Interval)`) rather than AoG's typical two-column
`mapping(:x, :y)` pointlike form.

Tested against AlgebraOfGraphics v0.13.1; its
`aesthetic_mapping` API is not yet fully stable across releases, so pin or
check compatibility if you hit `MethodError`/`aesthetic_mapping` errors on a
different version.

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
- Requires Julia 1.10+.
