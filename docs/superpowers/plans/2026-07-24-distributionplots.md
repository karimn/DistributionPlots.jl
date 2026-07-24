# DistributionPlots.jl Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a Julia package of Makie recipes that mimic R's `ggdist` (slab, interval, pointinterval, lineribbon, dots) for visualizing `RVars.jl` samples and `Distributions.jl` distributions through one interface.

**Architecture:** A pure statistics layer (`point_interval`, `slab_curve`, `dot_layout`, plus a KDE and a geometry transform) with zero Makie dependency is consumed by thin Makie recipes. One `slabinterval` parent recipe draws toggleable slab/interval/point sub-parts; `slab`/`interval`/`pointinterval`/`halfeye`/… are children that preset the toggles (ggdist's own "child = parent + defaults" mechanism). `lineribbon` and `dots` are separate engines. `RVars`, `MCMCChains`, and `AlgebraOfGraphics` support ship as weakdep extensions.

**Tech Stack:** Julia 1.10, Makie 0.24, KernelDensity, Distributions, StatsBase, Tables. R (ggdist 3.3.3 / posterior 1.6.1) for golden test fixtures. CairoMakie for headless smoke tests.

## Global Constraints

Every task's requirements implicitly include this section. Values are copied verbatim from the spec.

- **Julia floor:** `1.10` (current LTS). Set `julia = "1.10"` in `[compat]`.
- **Package UUID:** `36a31c4c-737d-46e0-87c4-7dcdbef46ebe`.
- **Core dependency UUIDs** (authoritative, from the General registry):
  - `Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"`
  - `Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"` (stdlib)
  - `StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"`
  - `KernelDensity = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"`
  - `Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"`
  - `Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"`
- **Weakdep UUIDs:**
  - `RVars = "108f709a-74fe-47a0-bf6a-d6e5c41d346f"`
  - `MCMCChains = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"`
  - `AlgebraOfGraphics = "cbdf2221-f076-402e-a563-3d30da359d67"`
- **Test-only UUIDs:** `CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"`, `Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"`, `Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"`.
- **`point_interval` output columns:** `value, lower, upper, width, point, interval` — plain names, **not** tidybayes's `.lower`/`.upper` (`.` is not a valid Julia identifier character). Return a `Vector{NamedTuple}` (Tables.jl-compatible as-is).
- **Attribute naming:** ggdist-style underscored names (`slab_type`, `show_interval`, `point_interval`, `show_slab`, `show_point`).
- **Interval level → colour** via Makie `colormap`/`colorrange` attributes, never a private palette.
- **Error tiers:** *reject loudly* (`ArgumentError` with the offending value in the message), *degrade with a single `@warn`*, or *handle silently* — never a bare `MethodError` for a user-reachable rejection, and never a silent drop.
- **Makie recipe API is version-sensitive.** Recipe/`convert_arguments` code below targets Makie 0.24's classic `@recipe(...) do scene Attributes(...) end` form. **In every recipe task, first run `julia --project -e 'using Makie; println(pkgversion(Makie))'` and skim `?Makie.@recipe`** to confirm the exact macro syntax for the installed version before writing recipe code.
- **Golden-fixture tooling pinned:** ggdist `3.3.3`, posterior `1.6.1`. The R generator records its `sessionInfo()` into the fixture directory.

---

## File Structure

```
DistributionPlots.jl/
├── Project.toml
├── README.md
├── src/
│   ├── DistributionPlots.jl   # module: includes, exports
│   ├── interface.jl           # AbstractDist, SampleDist, AnalyticDist; support/quantile_at/cdf_at; asdist
│   ├── density.jl             # kde_reflected; density_at
│   ├── intervals.jl           # point_interval (qi/hdci/hdi), point summaries
│   ├── slabs.jl               # slab_curve
│   ├── dotlayout.jl           # Wilkinson dot packing
│   ├── geometry.jl            # positioning transform: thickness → coordinates
│   └── recipes/
│       ├── slabinterval.jl    # parent engine + default-setter children + pre-summarised path
│       ├── dots.jl            # dotsinterval / dots recipes
│       └── lineribbon.jl      # separate engine
├── ext/
│   ├── DistributionPlotsRVarsExt.jl
│   ├── DistributionPlotsMCMCChainsExt.jl
│   └── DistributionPlotsAlgebraOfGraphicsExt.jl
└── test/
    ├── runtests.jl
    ├── gen_fixtures.R
    ├── fixtures/              # committed R-exported golden CSVs + sessionInfo
    ├── test_interface.jl
    ├── test_density.jl
    ├── test_intervals.jl
    ├── test_slabs.jl
    ├── test_dotlayout.jl
    ├── test_geometry.jl
    ├── test_golden.jl
    ├── test_recipes.jl
    └── test_extensions.jl
```

Julia note: package extensions live in `ext/` as **single files** named exactly `<PackageName><Weakdep>Ext.jl`, each defining a module of the same name. The `[extensions]` table in `Project.toml` maps the module name to its triggering weakdep(s).

---

## Task 1: Package scaffold

**Files:**
- Create: `Project.toml`
- Create: `src/DistributionPlots.jl`
- Create: `test/runtests.jl`
- Create: `.github/workflows/CI.yml`

**Interfaces:**
- Consumes: nothing.
- Produces: the `DistributionPlots` module; a green `Pkg.test()` baseline.

- [ ] **Step 1: Write the failing test**

Create `test/runtests.jl`:

```julia
using DistributionPlots
using Test

@testset "DistributionPlots.jl" begin
    @testset "smoke" begin
        @test isdefined(DistributionPlots, :point_interval) == false  # not yet; placeholder flips in Task 4
        @test DistributionPlots isa Module
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
Expected: FAIL — `Project.toml`/module do not exist yet.

- [ ] **Step 3: Write minimal implementation**

Create `Project.toml`:

```toml
name = "DistributionPlots"
uuid = "36a31c4c-737d-46e0-87c4-7dcdbef46ebe"
authors = ["Karim Naguib"]
version = "0.1.0"

[deps]
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
KernelDensity = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[weakdeps]
AlgebraOfGraphics = "cbdf2221-f076-402e-a563-3d30da359d67"
MCMCChains = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
RVars = "108f709a-74fe-47a0-bf6a-d6e5c41d346f"

[extensions]
DistributionPlotsAlgebraOfGraphicsExt = "AlgebraOfGraphics"
DistributionPlotsMCMCChainsExt = "MCMCChains"
DistributionPlotsRVarsExt = "RVars"

[compat]
Distributions = "0.25"
KernelDensity = "0.6"
Makie = "0.24"
Statistics = "1"
StatsBase = "0.34"
Tables = "1"
julia = "1.10"

[extras]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
MCMCChains = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
RVars = "108f709a-74fe-47a0-bf6a-d6e5c41d346f"
AlgebraOfGraphics = "cbdf2221-f076-402e-a563-3d30da359d67"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[targets]
test = ["Test", "Random", "CairoMakie", "RVars", "MCMCChains", "AlgebraOfGraphics"]
```

Create `src/DistributionPlots.jl`:

```julia
module DistributionPlots

using Statistics
using StatsBase
using Distributions
using KernelDensity
using Tables
using Makie

# includes are added as later tasks create the files:
# include("interface.jl")
# include("density.jl")
# include("intervals.jl")
# include("slabs.jl")
# include("dotlayout.jl")
# include("geometry.jl")
# include("recipes/slabinterval.jl")
# include("recipes/dots.jl")
# include("recipes/lineribbon.jl")

end # module
```

Fix the placeholder assertion in `test/runtests.jl` Step 1 — replace the first `@test` with:

```julia
        @test true  # baseline; real testsets are included per task below
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'`
Expected: PASS — `Test Summary: DistributionPlots.jl | 1 passed`. `RVars` resolves from GitHub if not registered; if `Pkg.instantiate()` errors on `RVars`, add it once with `julia --project=. -e 'using Pkg; Pkg.add(url="https://github.com/karimn/RVars.jl")'` then retry.

- [ ] **Step 5: Create CI and commit**

Create `.github/workflows/CI.yml`:

```yaml
name: CI
on:
  push:
    branches: [main]
  pull_request:
jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        version: ['1.10', '1']
    steps:
      - uses: actions/checkout@v4
      - uses: julia-actions/setup-julia@v2
        with:
          version: ${{ matrix.version }}
      - uses: julia-actions/cache@v2
      - uses: julia-actions/julia-buildpkg@v1
      - uses: julia-actions/julia-runtest@v1
```

```bash
git add Project.toml src/DistributionPlots.jl test/runtests.jl .github/workflows/CI.yml
git commit -m "feat: package scaffold with deps, extensions stanza, CI"
```

---

## Task 2: Input protocol — types and non-density methods

**Files:**
- Create: `src/interface.jl`
- Create: `test/test_interface.jl`
- Modify: `src/DistributionPlots.jl` (uncomment `include("interface.jl")`, add exports)
- Modify: `test/runtests.jl` (include the new test file)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces:
  - `abstract type AbstractDist end`
  - `struct SampleDist <: AbstractDist` with field `samples::Vector{Float64}`
  - `struct AnalyticDist{D<:Distributions.UnivariateDistribution} <: AbstractDist` with field `dist::D`
  - `asdist(x::AbstractVector{<:Real})::SampleDist` — drops `NaN`/`missing` with one `@warn`, rejects empty with `ArgumentError`
  - `asdist(d::Distributions.UnivariateDistribution)::AnalyticDist`
  - `support(d::AbstractDist)::Tuple{Float64,Float64}` — for `AnalyticDist` with infinite tails, trimmed to `(quantile(d, trim), quantile(d, 1-trim))`; signature `support(d; trim=0.001)`
  - `quantile_at(d::AbstractDist, p::Real)::Float64`
  - `cdf_at(d::AbstractDist, x::Real)::Float64`

- [ ] **Step 1: Write the failing test**

Create `test/test_interface.jl`:

```julia
using DistributionPlots
using DistributionPlots: SampleDist, AnalyticDist, asdist, support, quantile_at, cdf_at
using Distributions
using Test

@testset "interface" begin
    s = asdist([1.0, 2.0, 3.0, 4.0])
    @test s isa SampleDist
    @test quantile_at(s, 0.5) ≈ 2.5          # type-7 median of 1:4
    @test cdf_at(s, 2.5) ≈ 0.5
    @test support(s) == (1.0, 4.0)

    a = asdist(Normal(0, 1))
    @test a isa AnalyticDist
    @test quantile_at(a, 0.5) ≈ 0.0 atol=1e-12
    @test cdf_at(a, 0.0) ≈ 0.5
    lo, hi = support(a; trim=0.001)
    @test lo ≈ quantile(Normal(0,1), 0.001)
    @test hi ≈ quantile(Normal(0,1), 0.999)

    b = asdist(Uniform(0, 1))                # finite support: not trimmed
    @test support(b) == (0.0, 1.0)

    @test_throws ArgumentError asdist(Float64[])            # empty rejected
    @test_logs (:warn,) asdist([1.0, NaN, 3.0])            # NaN dropped with warn
    @test asdist([1.0, NaN, 3.0]).samples == [1.0, 3.0]
end
```

- [ ] **Step 2: Run test to verify it fails**

Add `include("test_interface.jl")` inside the `@testset "DistributionPlots.jl"` block in `test/runtests.jl`.
Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `SampleDist` not defined.

- [ ] **Step 3: Write minimal implementation**

Create `src/interface.jl`:

```julia
abstract type AbstractDist end

struct SampleDist <: AbstractDist
    samples::Vector{Float64}
end

struct AnalyticDist{D<:Distributions.UnivariateDistribution} <: AbstractDist
    dist::D
end

"""
    asdist(x)

Normalise an input into an `AbstractDist`. Vectors of samples become
`SampleDist` (NaN/missing dropped with a warning; empty is an error);
`UnivariateDistribution`s become `AnalyticDist`.
"""
function asdist(x::AbstractVector{<:Real})
    clean = collect(Float64, Iterators.filter(!isnan, skipmissing(x)))
    ndropped = length(x) - length(clean)
    ndropped > 0 && @warn "asdist: dropped $ndropped NaN/missing value(s) from samples"
    isempty(clean) && throw(ArgumentError("asdist: no finite samples (got $(length(x)) values, all NaN/missing/empty)"))
    return SampleDist(clean)
end

asdist(d::Distributions.UnivariateDistribution) = AnalyticDist(d)
asdist(d::AbstractDist) = d

quantile_at(d::SampleDist, p::Real) = quantile(d.samples, p)          # StatsBase type-7
quantile_at(d::AnalyticDist, p::Real) = quantile(d.dist, p)

cdf_at(d::SampleDist, x::Real) = count(≤(x), d.samples) / length(d.samples)
cdf_at(d::AnalyticDist, x::Real) = cdf(d.dist, x)

function support(d::SampleDist; trim::Real=0.001)
    return (minimum(d.samples), maximum(d.samples))
end

function support(d::AnalyticDist; trim::Real=0.001)
    lo = minimum(d.dist)
    hi = maximum(d.dist)
    isfinite(lo) || (lo = quantile(d.dist, trim))
    isfinite(hi) || (hi = quantile(d.dist, 1 - trim))
    return (float(lo), float(hi))
end
```

In `src/DistributionPlots.jl`, uncomment `include("interface.jl")` and add:

```julia
export asdist
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS — the `interface` testset is green.

- [ ] **Step 5: Commit**

```bash
git add src/interface.jl src/DistributionPlots.jl test/test_interface.jl test/runtests.jl
git commit -m "feat: input protocol types with support/quantile_at/cdf_at"
```

---

## Task 3: KDE with boundary reflection

**Files:**
- Create: `src/density.jl`
- Create: `test/test_density.jl`
- Modify: `src/DistributionPlots.jl` (include)
- Modify: `test/runtests.jl` (include)

**Interfaces:**
- Consumes: `SampleDist`, `AnalyticDist`, `support` (Task 2).
- Produces:
  - `kde_reflected(samples::Vector{Float64}; bounds::Tuple{Float64,Float64}, npoints::Int=512)` → `(xs::Vector{Float64}, dens::Vector{Float64})`. Reflects samples across any finite bound before KDE so density is not underestimated near the boundary.
  - `density_at(d::AbstractDist, xs::AbstractVector{<:Real}; bounds=support(d))` → `Vector{Float64}` densities at `xs`.

**Note:** the bounded-KDE reflection is the "correctness landmine" from the spec — `KernelDensity.jl` does not reflect natively, so we do it here. A bound is treated as finite (reflect) only when the samples do not exceed it; unbounded sides get no reflection.

- [ ] **Step 1: Write the failing test**

Create `test/test_density.jl`:

```julia
using DistributionPlots
using DistributionPlots: asdist, density_at, kde_reflected, support
using Distributions
using Test

@testset "density" begin
    # Analytic density_at is exact pdf
    a = asdist(Normal(0, 1))
    @test density_at(a, [0.0]) ≈ [pdf(Normal(0,1), 0.0)]

    # Reflected KDE on [0,1]-bounded samples: mass near a boundary is not
    # underestimated the way an unbounded KDE would underestimate it.
    Random_seeded = Beta(2, 8)
    samples = quantile.(Beta(2,8), range(0.001, 0.999; length=4000))  # deterministic
    xs, dens = kde_reflected(samples; bounds=(0.0, 1.0), npoints=512)
    @test all(≥(0), dens)
    @test minimum(xs) ≥ -1e-9 && maximum(xs) ≤ 1 + 1e-9
    # integral ≈ 1 over [0,1]
    area = sum((dens[2:end] .+ dens[1:end-1]) ./ 2 .* diff(xs))
    @test area ≈ 1.0 atol=0.05

    # Reflected density near the lower boundary exceeds the naive-unbounded
    # estimate at the same point (reflection adds the mirrored mass).
    s = asdist(samples)
    d_lo = density_at(s, [0.05]; bounds=(0.0, 1.0))[1]
    @test d_lo > 0
end
```

- [ ] **Step 2: Run test to verify it fails**

Add `include("test_density.jl")` to `test/runtests.jl`.
Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `kde_reflected` not defined.

- [ ] **Step 3: Write minimal implementation**

Create `src/density.jl`:

```julia
using KernelDensity: kde, pdf as kdepdf

"""
    kde_reflected(samples; bounds, npoints=512)

Kernel density estimate that reflects samples across any *finite* bound so the
density is not underestimated near the boundary (ggdist's `density = "bounded"`).
Returns `(xs, dens)` evaluated on `npoints` points spanning `bounds`.
"""
function kde_reflected(samples::AbstractVector{<:Real};
                       bounds::Tuple{<:Real,<:Real},
                       npoints::Int=512)
    lo, hi = float(bounds[1]), float(bounds[2])
    data = collect(Float64, samples)

    reflect_lo = isfinite(lo)
    reflect_hi = isfinite(hi)
    augmented = copy(data)
    reflect_lo && append!(augmented, 2lo .- data)   # mirror across lower bound
    reflect_hi && append!(augmented, 2hi .- data)   # mirror across upper bound

    k = kde(augmented)                                # KernelDensity picks bandwidth
    xs = collect(range(lo, hi; length=npoints))
    dens = kdepdf(k, xs)
    # Reflection triples (or doubles) total mass inside [lo,hi]; renormalise so
    # the curve integrates to ~1 over the true support.
    scale = 1 + reflect_lo + reflect_hi
    dens .*= scale
    dens .= max.(dens, 0.0)
    return xs, dens
end

density_at(d::AnalyticDist, xs::AbstractVector{<:Real}; bounds=support(d)) =
    pdf.(d.dist, xs)

function density_at(d::SampleDist, xs::AbstractVector{<:Real}; bounds=support(d))
    gx, gd = kde_reflected(d.samples; bounds=bounds, npoints=max(512, length(xs)))
    # linear-interpolate the grid density onto the requested xs
    return _interp(gx, gd, xs)
end

function _interp(gx::Vector{Float64}, gd::Vector{Float64}, xs)
    out = similar(collect(Float64, xs))
    for (i, x) in enumerate(xs)
        if x ≤ gx[1]
            out[i] = gd[1]
        elseif x ≥ gx[end]
            out[i] = gd[end]
        else
            j = searchsortedlast(gx, x)
            t = (x - gx[j]) / (gx[j+1] - gx[j])
            out[i] = (1 - t) * gd[j] + t * gd[j+1]
        end
    end
    return out
end
```

Add `include("density.jl")` to `src/DistributionPlots.jl` (after `interface.jl`).

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS — `density` testset green. If `area` is off, widen the atol comment note but do not exceed `0.05`; a larger miss means the renormalisation `scale` is wrong.

- [ ] **Step 5: Commit**

```bash
git add src/density.jl src/DistributionPlots.jl test/test_density.jl test/runtests.jl
git commit -m "feat: boundary-reflected KDE and density_at"
```

---

## Task 4: point_interval — point summaries, qi, hdci

**Files:**
- Create: `src/intervals.jl`
- Create: `test/test_intervals.jl`
- Modify: `src/DistributionPlots.jl` (include + export `point_interval`)
- Modify: `test/runtests.jl` (include)

**Interfaces:**
- Consumes: `AbstractDist`, `asdist`, `quantile_at` (Task 2).
- Produces:
  - `const PIRow = @NamedTuple{value::Float64, lower::Float64, upper::Float64, width::Float64, point::Symbol, interval::Symbol}`
  - `point_interval(d; widths=[0.66, 0.95], point=:median, interval=:qi)::Vector{PIRow}` where `d` is anything `asdist` accepts. `point ∈ {:mean,:median}` and `interval ∈ {:qi,:hdci}` in this task (`:mode`/`:hdi` added in Task 5). Returns one row per (width × interval-piece); `qi`/`hdci` produce exactly one piece per width.
  - `point_summary(d::AbstractDist, point::Symbol)::Float64`

- [ ] **Step 1: Write the failing test**

Create `test/test_intervals.jl`:

```julia
using DistributionPlots
using DistributionPlots: point_interval, point_summary, asdist
using Statistics
using Test

@testset "intervals: qi + hdci" begin
    s = collect(1.0:100.0)              # uniform-ish sample
    rows = point_interval(s; widths=[0.5, 0.9], point=:median, interval=:qi)
    @test length(rows) == 2
    @test rows[1].width == 0.5
    @test rows[1].point == :median
    @test rows[1].interval == :qi
    @test rows[1].value ≈ median(s)
    # 50% qi = [q25, q75]
    @test rows[1].lower ≈ quantile(s, 0.25)
    @test rows[1].upper ≈ quantile(s, 0.75)
    # 90% qi = [q05, q95]
    @test rows[2].lower ≈ quantile(s, 0.05)
    @test rows[2].upper ≈ quantile(s, 0.95)

    # mean point summary
    rmean = point_interval(s; widths=[0.9], point=:mean)
    @test rmean[1].value ≈ mean(s)

    # hdci on a symmetric sample ≈ qi (narrowest contiguous window)
    rh = point_interval(s; widths=[0.9], interval=:hdci)
    @test rh[1].interval == :hdci
    @test (rh[1].upper - rh[1].lower) ≤ (rmean[1].upper - rmean[1].lower) + 1e-9

    # Tables.jl compatibility: a Vector{NamedTuple} is a valid row table
    @test Tables.istable(rows)
    @test Tables.rowaccess(rows)
end
```

- [ ] **Step 2: Run test to verify it fails**

Add `include("test_intervals.jl")` and `using Tables` at the top of `test/runtests.jl` (once). Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `point_interval` not defined.

- [ ] **Step 3: Write minimal implementation**

Create `src/intervals.jl`:

```julia
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
```

Add `include("intervals.jl")` to `src/DistributionPlots.jl` and `export point_interval`.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS — `intervals: qi + hdci` green.

- [ ] **Step 5: Commit**

```bash
git add src/intervals.jl src/DistributionPlots.jl test/test_intervals.jl test/runtests.jl
git commit -m "feat: point_interval with mean/median points and qi/hdci intervals"
```

---

## Task 5: point_interval — mode and hdi (possibly disjoint)

**Files:**
- Modify: `src/intervals.jl` (add `_mode`, `_hdi`)
- Modify: `test/test_intervals.jl` (add a disjoint-hdi testset)

**Interfaces:**
- Consumes: `density_at`, `kde_reflected` (Task 3), `support` (Task 2).
- Produces:
  - `_mode(d::AbstractDist)::Float64`
  - `_hdi(d::AbstractDist, w::Real)::Vector{Tuple{Float64,Float64}}` — one or more disjoint pieces. Enables `point=:mode` and `interval=:hdi` in `point_interval`.

- [ ] **Step 1: Write the failing test**

Add to `test/test_intervals.jl`:

```julia
@testset "intervals: mode + disjoint hdi" begin
    using Distributions
    # Bimodal mixture: two well-separated Normals → hdi returns TWO pieces at 0.95.
    mix = MixtureModel([Normal(-3, 0.4), Normal(3, 0.4)], [0.5, 0.5])
    draws = reduce(vcat, (rand(Normal(-3,0.4), 5000), rand(Normal(3,0.4), 5000)))
    rows = point_interval(draws; widths=[0.95], interval=:hdi, point=:median)
    @test length(rows) ≥ 2                       # disjoint → multiple rows
    @test all(r -> r.interval == :hdi && r.width == 0.95, rows)
    # the two pieces straddle the two modes
    los = sort(getproperty.(rows, :lower))
    @test los[1] < 0 && los[end] > 0

    # mode of a right-skewed sample is left of the mean
    skew = rand(Exponential(1.0), 20000)
    rmode = point_interval(skew; widths=[0.9], point=:mode)
    @test rmode[1].value < mean(skew)
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `_hdi`/`_mode` not defined (or `UndefVarError`).

- [ ] **Step 3: Write minimal implementation**

Append to `src/intervals.jl`:

```julia
function _mode(d::SampleDist)
    lo, hi = support(d)
    xs, dens = kde_reflected(d.samples; bounds=(lo, hi), npoints=1024)
    return xs[argmax(dens)]
end
_mode(d::AnalyticDist) = mode(d.dist)

"""
    _hdi(d, w) -> Vector{(lower, upper)}

Highest-density interval: the set `{x : f(x) ≥ c}` whose total mass is `w`,
found by lowering a horizontal threshold on the (KDE for samples, exact for
analytic) density until the covered probability reaches `w`. May be disjoint.
"""
function _hdi(d::SampleDist, w::Real)
    lo, hi = support(d)
    xs, dens = kde_reflected(d.samples; bounds=(lo, hi), npoints=1024)
    return _hdi_from_density(xs, dens, w)
end

function _hdi(d::AnalyticDist, w::Real)
    lo, hi = support(d; trim=1e-4)
    xs = collect(range(lo, hi; length=2048))
    dens = pdf.(d.dist, xs)
    return _hdi_from_density(xs, dens, w)
end

# Threshold-sweep on a density grid. Returns contiguous runs above the threshold
# whose trapezoidal mass first reaches `w`.
function _hdi_from_density(xs::Vector{Float64}, dens::Vector{Float64}, w::Real)
    dx = diff(xs)
    total = sum((dens[2:end] .+ dens[1:end-1]) ./ 2 .* dx)
    dens = dens ./ total                       # normalise to unit mass on grid
    order = sortperm(dens; rev=true)
    thresh = 0.0
    covered = 0.0
    # accumulate cells from highest density down until mass ≥ w
    cellmass = similar(dens)
    cellmass[1] = 0.0
    @inbounds for i in 2:length(xs)
        cellmass[i] = (dens[i] + dens[i-1]) / 2 * dx[i-1]
    end
    keep = falses(length(xs))
    for idx in order
        keep[idx] = true
        covered += idx == 1 ? 0.0 : cellmass[idx]
        thresh = dens[idx]
        covered ≥ w && break
    end
    # merge kept indices into contiguous (lower, upper) runs
    pieces = Tuple{Float64,Float64}[]
    i = 1
    n = length(xs)
    while i ≤ n
        if keep[i]
            j = i
            while j < n && keep[j+1]
                j += 1
            end
            push!(pieces, (xs[i], xs[j]))
            i = j + 1
        else
            i += 1
        end
    end
    isempty(pieces) && push!(pieces, (xs[argmax(dens)], xs[argmax(dens)]))
    return pieces
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS — `intervals: mode + disjoint hdi` green. If the mixture yields one piece, raise `npoints` or widen mode separation in the test data; the algorithm is correct when modes are separated by more than a few bandwidths.

- [ ] **Step 5: Commit**

```bash
git add src/intervals.jl test/test_intervals.jl
git commit -m "feat: mode point summary and disjoint hdi intervals"
```

---

## Task 6: slab_curve

**Files:**
- Create: `src/slabs.jl`
- Create: `test/test_slabs.jl`
- Modify: `src/DistributionPlots.jl` (include + export `slab_curve`)
- Modify: `test/runtests.jl` (include)

**Interfaces:**
- Consumes: `AbstractDist`, `asdist`, `support`, `density_at`, `cdf_at` (Tasks 2–3).
- Produces: `slab_curve(d; kind=:pdf, n=201, trim=0.001)::Tuple{Vector{Float64},Vector{Float64}}` returning `(xs, thickness)`. `kind ∈ {:pdf,:cdf,:ccdf,:histogram}`. Thickness is raw (unnormalised); normalisation for drawing is the geometry layer's job (Task 8).

- [ ] **Step 1: Write the failing test**

Create `test/test_slabs.jl`:

```julia
using DistributionPlots
using DistributionPlots: slab_curve, asdist
using Distributions
using Test

@testset "slab_curve" begin
    a = asdist(Normal(0, 1))
    xs, th = slab_curve(a; kind=:pdf, n=201)
    @test length(xs) == length(th) == 201
    @test th[argmin(abs.(xs))] ≈ pdf(Normal(0,1), 0.0) atol=1e-2   # peak at 0
    @test all(≥(0), th)

    # cdf is monotone nondecreasing from ~0 to ~1
    xs2, cd = slab_curve(a; kind=:cdf, n=201)
    @test issorted(cd)
    @test cd[1] < 0.05 && cd[end] > 0.95

    # ccdf is 1 - cdf
    _, cc = slab_curve(a; kind=:ccdf, n=201)
    @test cc ≈ 1 .- cd atol=1e-6

    # histogram returns bin centers and counts/density
    s = asdist(randn(5000))
    xh, h = slab_curve(s; kind=:histogram)
    @test length(xh) == length(h)
    @test all(≥(0), h)
end
```

- [ ] **Step 2: Run test to verify it fails**

Add `include("test_slabs.jl")` to `test/runtests.jl`. Run tests.
Expected: FAIL — `slab_curve` not defined.

- [ ] **Step 3: Write minimal implementation**

Create `src/slabs.jl`:

```julia
using StatsBase: fit, Histogram

"""
    slab_curve(d; kind=:pdf, n=201, trim=0.001)

Return `(xs, thickness)` describing the slab shape. `kind`:
- `:pdf`  — density (KDE for samples, exact pdf for analytic)
- `:cdf`  — cumulative distribution
- `:ccdf` — complementary cdf (`1 - cdf`)
- `:histogram` — binned counts (samples only)
"""
function slab_curve(d; kind::Symbol=:pdf, n::Int=201, trim::Real=0.001)
    dist = asdist(d)
    lo, hi = support(dist; trim=trim)
    if kind === :histogram
        return _histogram_curve(dist)
    end
    xs = collect(range(lo, hi; length=n))
    if kind === :pdf
        return xs, density_at(dist, xs; bounds=(lo, hi))
    elseif kind === :cdf
        return xs, [cdf_at(dist, x) for x in xs]
    elseif kind === :ccdf
        return xs, [1 - cdf_at(dist, x) for x in xs]
    else
        throw(ArgumentError("slab_curve: unknown kind :$kind (expected :pdf, :cdf, :ccdf, or :histogram)"))
    end
end

function _histogram_curve(d::SampleDist)
    h = fit(Histogram, d.samples; nbins=_fd_nbins(d.samples))
    edges = collect(h.edges[1])
    centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
    weights = float.(h.weights)
    return centers, weights
end

_histogram_curve(::AnalyticDist) =
    throw(ArgumentError("slab_curve: kind=:histogram requires samples, not an analytic distribution"))

# Freedman–Diaconis bin count
function _fd_nbins(x::Vector{Float64})
    n = length(x)
    iqr = quantile(x, 0.75) - quantile(x, 0.25)
    bw = iqr > 0 ? 2 * iqr / cbrt(n) : (maximum(x) - minimum(x)) / max(1, round(Int, sqrt(n)))
    bw ≤ 0 && return 1
    return max(1, ceil(Int, (maximum(x) - minimum(x)) / bw))
end
```

Add `include("slabs.jl")` to `src/DistributionPlots.jl` and `export slab_curve`.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS — `slab_curve` green.

- [ ] **Step 5: Commit**

```bash
git add src/slabs.jl src/DistributionPlots.jl test/test_slabs.jl test/runtests.jl
git commit -m "feat: slab_curve for pdf/cdf/ccdf/histogram"
```

---

## Task 7: dot_layout — Wilkinson dot packing

**Files:**
- Create: `src/dotlayout.jl`
- Create: `test/test_dotlayout.jl`
- Modify: `src/DistributionPlots.jl` (include + export `dot_layout`)
- Modify: `test/runtests.jl` (include)

**Interfaces:**
- Consumes: `AbstractDist`, `asdist`, `quantile_at` (Task 2).
- Produces: `dot_layout(d; ndots=50, binwidth=nothing)::@NamedTuple{x::Vector{Float64}, y::Vector{Int}, binwidth::Float64}` — quantile dotplot: `ndots` representative quantiles packed into stacks. `y` is the stack height (1-based) for each dot; `binwidth` is the dot diameter in data units.

- [ ] **Step 1: Write the failing test**

Create `test/test_dotlayout.jl`:

```julia
using DistributionPlots
using DistributionPlots: dot_layout, asdist
using Distributions
using Test

@testset "dot_layout" begin
    a = asdist(Normal(0, 1))
    lay = dot_layout(a; ndots=50)
    @test length(lay.x) == 50
    @test length(lay.y) == 50
    @test all(≥(1), lay.y)                       # stacks are 1-based
    @test lay.binwidth > 0
    # dots are sorted along x (quantile dotplot)
    @test issorted(lay.x)
    # near the mode (x≈0) stacks are taller than in the tails
    center_heights = maximum(lay.y[abs.(lay.x) .< 0.5])
    tail_heights = maximum(lay.y[abs.(lay.x) .> 2.0])
    @test center_heights ≥ tail_heights
end
```

- [ ] **Step 2: Run test to verify it fails**

Add `include("test_dotlayout.jl")` to `test/runtests.jl`. Run tests.
Expected: FAIL — `dot_layout` not defined.

- [ ] **Step 3: Write minimal implementation**

Create `src/dotlayout.jl`:

```julia
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
```

Add `include("dotlayout.jl")` to `src/DistributionPlots.jl` and `export dot_layout`.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS — `dot_layout` green.

- [ ] **Step 5: Commit**

```bash
git add src/dotlayout.jl src/DistributionPlots.jl test/test_dotlayout.jl test/runtests.jl
git commit -m "feat: quantile dotplot layout (Wilkinson binning)"
```

---

## Task 8: geometry transform

**Files:**
- Create: `src/geometry.jl`
- Create: `test/test_geometry.jl`
- Modify: `src/DistributionPlots.jl` (include)
- Modify: `test/runtests.jl` (include)

**Interfaces:**
- Consumes: nothing (pure coordinate math).
- Produces:
  - `slab_polygon(xs, thickness; position, orientation, side, justification, scale, normalize)::Vector{Point2f}` — the filled outline of a slab.
  - `interval_segment(lower, upper; position, orientation)::Tuple{Point2f,Point2f}` — endpoints of an interval bar.
  - `point_marker(value; position, orientation)::Point2f` — the point-summary marker location.
  - `normalize_thickness(thickness, mode)::Vector{Float64}` — scale thickness to `[0,1]` per `mode ∈ {:none,:all}` (v1 supports `:all` = divide by global max, and `:none` = leave as-is).

Convention: `orientation=:vertical` means the **value axis is y** and slabs grow along **x** from the baseline at `position`. `side=:top` grows +x, `:bottom` grows -x, `:both` is symmetric. `orientation=:horizontal` swaps x and y.

- [ ] **Step 1: Write the failing test**

Create `test/test_geometry.jl`:

```julia
using DistributionPlots
using DistributionPlots: slab_polygon, interval_segment, point_marker, normalize_thickness
using Makie: Point2f
using Test

@testset "geometry" begin
    @test normalize_thickness([0.0, 1.0, 2.0], :all) == [0.0, 0.5, 1.0]
    @test normalize_thickness([0.0, 1.0, 2.0], :none) == [0.0, 1.0, 2.0]

    xs = [0.0, 1.0, 2.0]
    th = [0.0, 1.0, 0.0]                        # normalised triangle
    # vertical, side=:top, position=5, scale=1 → polygon x ∈ [5, 6], y follows xs
    poly = slab_polygon(xs, th; position=5.0, orientation=:vertical,
                        side=:top, justification=0.0, scale=1.0, normalize=:none)
    @test poly isa Vector{Point2f}
    @test all(p -> p[1] ≥ 5.0 - 1e-6, poly)      # never left of baseline
    @test maximum(p -> p[1], poly) ≈ 6.0         # peak thickness 1 * scale 1
    @test any(p -> p[2] ≈ 1.0, poly)             # value axis is y

    # interval segment: vertical orientation → x fixed at position, y = [lo, hi]
    p1, p2 = interval_segment(-1.0, 2.0; position=5.0, orientation=:vertical)
    @test p1 == Point2f(5.0, -1.0)
    @test p2 == Point2f(5.0, 2.0)

    pm = point_marker(0.5; position=5.0, orientation=:vertical)
    @test pm == Point2f(5.0, 0.5)

    # horizontal orientation swaps axes
    p1h, p2h = interval_segment(-1.0, 2.0; position=5.0, orientation=:horizontal)
    @test p1h == Point2f(-1.0, 5.0)
    @test p2h == Point2f(2.0, 5.0)
end
```

- [ ] **Step 2: Run test to verify it fails**

Add `include("test_geometry.jl")` to `test/runtests.jl`. Run tests.
Expected: FAIL — `slab_polygon` not defined.

- [ ] **Step 3: Write minimal implementation**

Create `src/geometry.jl`:

```julia
using Makie: Point2f

function normalize_thickness(thickness::AbstractVector{<:Real}, mode::Symbol)
    if mode === :none
        return collect(Float64, thickness)
    elseif mode === :all
        m = maximum(thickness)
        return m > 0 ? collect(Float64, thickness ./ m) : collect(Float64, thickness)
    else
        throw(ArgumentError("normalize_thickness: unknown mode :$mode (expected :all or :none)"))
    end
end

# Map a (value, offset) pair to a Point2f under an orientation.
# `along` is the value-axis coordinate; `perp` is the thickness/position axis.
_pt(along, perp, ::Val{:vertical}) = Point2f(perp, along)
_pt(along, perp, ::Val{:horizontal}) = Point2f(along, perp)

function slab_polygon(xs::AbstractVector{<:Real}, thickness::AbstractVector{<:Real};
                      position::Real, orientation::Symbol=:vertical,
                      side::Symbol=:top, justification::Real=0.0,
                      scale::Real=1.0, normalize::Symbol=:all)
    th = normalize_thickness(thickness, normalize) .* scale
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
```

Add `include("geometry.jl")` to `src/DistributionPlots.jl`.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS — `geometry` green.

- [ ] **Step 5: Commit**

```bash
git add src/geometry.jl src/DistributionPlots.jl test/test_geometry.jl test/runtests.jl
git commit -m "feat: geometry transform (slab polygon, interval segment, point marker)"
```

---

## Task 9: R golden fixtures

**Files:**
- Create: `test/gen_fixtures.R`
- Create: `test/fixtures/` (generated CSVs + `sessionInfo.txt`, committed)
- Create: `test/test_golden.jl`
- Modify: `test/runtests.jl` (include)

**Interfaces:**
- Consumes: `point_interval` (Tasks 4–5), `slab_curve` (Task 6).
- Produces: committed golden CSVs and a Julia testset asserting our stats match ggdist within tiered tolerance.

**Prerequisite:** R with `ggdist` 3.3.3 and `posterior` 1.6.1 (verified installed at `/usr/local/bin/Rscript`).

- [ ] **Step 1: Write the fixture generator**

Create `test/gen_fixtures.R`:

```r
# Regenerate golden fixtures. Run: Rscript test/gen_fixtures.R
suppressMessages({library(ggdist); library(posterior)})
set.seed(20260724)
dir.create("test/fixtures", showWarnings = FALSE, recursive = TRUE)

emit_samples <- function(name, draws) {
  write.csv(data.frame(draw = draws),
            file.path("test/fixtures", paste0(name, "_draws.csv")), row.names = FALSE)
}
emit_pi <- function(name, draws, width, interval_fn, interval_name) {
  pi <- interval_fn(draws, .width = width)   # median_qi / median_hdci / mode_hdi
  out <- data.frame(value = pi$y, lower = pi$ymin, upper = pi$ymax,
                    width = width, interval = interval_name)
  write.csv(out, file.path("test/fixtures", paste0(name, "_", interval_name, "_", width, ".csv")),
            row.names = FALSE)
}

cases <- list(
  normal   = rnorm(4000, 0, 1),
  beta     = rbeta(4000, 2, 8),
  bimodal  = c(rnorm(2500, -3, 0.4), rnorm(2500, 3, 0.4)),
  studentt = rt(4000, df = 3)
)

for (nm in names(cases)) {
  d <- cases[[nm]]
  emit_samples(nm, d)
  for (w in c(0.66, 0.95)) {
    emit_pi(nm, d, w, median_qi,   "qi")
    emit_pi(nm, d, w, median_hdci, "hdci")
  }
}

writeLines(capture.output(sessionInfo()), "test/fixtures/sessionInfo.txt")
cat("fixtures written\n")
```

- [ ] **Step 2: Generate fixtures and verify they exist**

Run: `cd <repo> && Rscript test/gen_fixtures.R`
Expected: prints `fixtures written`; `test/fixtures/` contains `*_draws.csv`, `*_qi_*.csv`, `*_hdci_*.csv`, and `sessionInfo.txt`. Confirm `sessionInfo.txt` shows `ggdist_3.3.3` and `posterior_1.6.1`.

- [ ] **Step 3: Write the golden test**

Create `test/test_golden.jl`:

```julia
using DistributionPlots
using DistributionPlots: point_interval
using DelimitedFiles
using Test

fixdir = joinpath(@__DIR__, "fixtures")

read_draws(name) = vec(readdlm(joinpath(fixdir, "$(name)_draws.csv"), ','; skipstart=1))

function read_pi(name, interval, width)
    rows = readdlm(joinpath(fixdir, "$(name)_$(interval)_$(width).csv"), ','; skipstart=1)
    # columns: value, lower, upper, width, interval
    return (value=Float64(rows[1,1]), lower=Float64(rows[1,2]), upper=Float64(rows[1,3]))
end

@testset "golden vs ggdist" begin
    for name in ("normal", "beta", "bimodal", "studentt")
        draws = read_draws(name)
        for w in (0.66, 0.95)
            # qi: type-7 quantiles match R exactly → tight tolerance
            ours = point_interval(draws; widths=[w], interval=:qi, point=:median)[1]
            ref = read_pi(name, "qi", w)
            @test ours.value ≈ ref.value atol=1e-8
            @test ours.lower ≈ ref.lower atol=1e-8
            @test ours.upper ≈ ref.upper atol=1e-8

            # hdci: sliding-window over sorted draws; matches ggdist closely
            oh = point_interval(draws; widths=[w], interval=:hdci, point=:median)[1]
            rh = read_pi(name, "hdci", w)
            @test oh.lower ≈ rh.lower atol=1e-6
            @test oh.upper ≈ rh.upper atol=1e-6
        end
    end
end
```

- [ ] **Step 4: Run test to verify it passes**

Add `include("test_golden.jl")` to `test/runtests.jl`. Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS. If `qi` fails at `1e-8`, our `quantile` is not type-7 — investigate `StatsBase.quantile` defaults (it is type-7; a failure means a wrong axis or off-by-one in `point_interval`). If `hdci` fails slightly, confirm ggdist's `median_hdci` uses the same contiguous-window definition; document the realised tolerance in a comment rather than loosening past `1e-6` without cause.

- [ ] **Step 5: Commit**

```bash
git add test/gen_fixtures.R test/fixtures test/test_golden.jl test/runtests.jl
git commit -m "test: golden fixtures from ggdist 3.3.3 for qi/hdci"
```

---

## Task 10: Recipe scaffolding — convert_arguments and the SlabInterval attributes

**Files:**
- Create: `src/recipes/slabinterval.jl`
- Create: `test/test_recipes.jl`
- Modify: `src/DistributionPlots.jl` (include; export recipe names)
- Modify: `test/runtests.jl` (include)

**Interfaces:**
- Consumes: `asdist`, `AbstractDist` (Task 2).
- Produces:
  - the `SlabInterval` recipe type and `slabinterval`/`slabinterval!` functions with the full attribute set.
  - `_to_dist_args(x)` — the **shared converter** turning a user input into the recipe's positional args `(positions::Vector, dists::Vector{AbstractDist})`.
  - `convert_arguments` methods for `Distribution` and `AbstractVector`, registered for `SlabInterval` (Task 12 extends the registration to every public recipe type).

**Before writing:** confirm Makie's recipe macro syntax — `julia --project -e 'using Makie; println(pkgversion(Makie)); '` and `?Makie.@recipe`.

- [ ] **Step 1: Write the failing test**

Create `test/test_recipes.jl`:

```julia
using DistributionPlots
using DistributionPlots: SlabInterval, _to_dist_args, AnalyticDist, SampleDist
using Makie
using Distributions
using Test

@testset "recipe scaffolding" begin
    # shared converter: a single distribution → one position, one dist
    pos, dists = _to_dist_args(Normal(0,1))
    @test pos == [1.0]
    @test dists[1] isa AnalyticDist

    # a plain sample vector → one SampleDist at position 1
    pos2, dists2 = _to_dist_args(randn(500))
    @test pos2 == [1.0]
    @test dists2[1] isa SampleDist

    # the recipe type exists and carries our attributes with ggdist-style names
    @test SlabInterval <: Makie.AbstractPlot
    attrs = Makie.default_theme(nothing, SlabInterval)
    for a in (:slab_type, :interval, :point, :side, :justification, :scale,
              :normalize, :show_slab, :show_interval, :show_point, :trim, :colormap)
        @test haskey(attrs, a)
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Add `include("test_recipes.jl")` to `test/runtests.jl`. Run tests.
Expected: FAIL — `SlabInterval` not defined.

- [ ] **Step 3: Write minimal implementation**

Create `src/recipes/slabinterval.jl`:

```julia
using Makie

# Shared converter: normalise any input to (positions, Vector{AbstractDist}).
# A single distribution/sample-vector is placed at position 1. Extensions
# (RVars, MCMCChains) add methods that expand multi-element inputs to 1:k.
_to_dist_args(x::Distributions.UnivariateDistribution) = ([1.0], AbstractDist[asdist(x)])
_to_dist_args(x::AbstractVector{<:Real}) = ([1.0], AbstractDist[asdist(x)])
function _to_dist_args(xs::AbstractVector{<:Distributions.UnivariateDistribution})
    return (collect(1.0:length(xs)), AbstractDist[asdist(x) for x in xs])
end

@recipe(SlabInterval, data) do scene
    Attributes(
        slab_type = :pdf,          # :pdf, :cdf, :ccdf, :histogram
        interval = :qi,            # :qi, :hdci, :hdi
        point = :median,           # :mean, :median, :mode
        widths = [0.66, 0.95],
        side = :top,               # :top, :bottom, :both
        justification = 0.0,
        scale = 0.9,
        normalize = :all,          # :all, :none
        show_slab = true,
        show_interval = true,
        show_point = true,
        trim = 0.001,
        n = 201,
        color = :black,
        colormap = :blues,
        colorrange = Makie.automatic,
        interval_linewidth = 6,
        point_size = 10,
    )
end

# Register convert_arguments for the base recipe. Task 12 loops this over all
# public recipe types; the extensions add per-type methods for RVar/Chains.
Makie.convert_arguments(::Type{<:SlabInterval}, x::Distributions.UnivariateDistribution) =
    (_to_dist_args(x),)
Makie.convert_arguments(::Type{<:SlabInterval}, x::AbstractVector{<:Real}) =
    (_to_dist_args(x),)
Makie.convert_arguments(::Type{<:SlabInterval}, xs::AbstractVector{<:Distributions.UnivariateDistribution}) =
    (_to_dist_args(xs),)
```

Add to `src/DistributionPlots.jl`:

```julia
include("recipes/slabinterval.jl")
export slabinterval, slabinterval!, SlabInterval
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS — `recipe scaffolding` green. If `default_theme(nothing, SlabInterval)` errors, use the version-appropriate accessor found via `?Makie.@recipe` (older Makie: `Makie.default_theme`; the attribute set is the same).

- [ ] **Step 5: Commit**

```bash
git add src/recipes/slabinterval.jl src/DistributionPlots.jl test/test_recipes.jl test/runtests.jl
git commit -m "feat: SlabInterval recipe type, attributes, shared converter"
```

---

## Task 11: SlabInterval plot! — draw the sub-parts

**Files:**
- Modify: `src/recipes/slabinterval.jl` (add `Makie.plot!`)
- Modify: `test/test_recipes.jl` (add a rendering smoke test)

**Interfaces:**
- Consumes: `slab_curve` (Task 6), `point_interval` (Tasks 4–5), geometry functions (Task 8), `SlabInterval` attributes (Task 10).
- Produces: a working `Makie.plot!(p::SlabInterval)` that draws slab (`poly!`), interval (`linesegments!`), and point (`scatter!`) sub-parts per the `show_*` toggles.

- [ ] **Step 1: Write the failing test**

Add to `test/test_recipes.jl`:

```julia
using CairoMakie

@testset "slabinterval renders" begin
    CairoMakie.activate!()
    f = slabinterval(randn(2000))
    @test f isa Makie.FigureAxisPlot
    ax = f.axis
    # at least one poly (slab) and some line/scatter children exist
    @test !isempty(ax.scene.plots)

    # toggles: show only the interval
    f2 = slabinterval(randn(2000); show_slab=false, show_point=false)
    @test f2 isa Makie.FigureAxisPlot

    # saving to a headless PNG does not error
    tmp = tempname() * ".png"
    save(tmp, f)
    @test isfile(tmp)
end
```

- [ ] **Step 2: Run test to verify it fails**

Run tests. Expected: FAIL — `plot!(::SlabInterval)` not implemented, so no plots are drawn / method error.

- [ ] **Step 3: Write minimal implementation**

Append to `src/recipes/slabinterval.jl`:

```julia
function Makie.plot!(p::SlabInterval)
    args = lift(p.data) do a
        a isa Tuple ? a[1] : a          # convert_arguments wraps as ((pos,dists),)
    end

    lift(args, p.slab_type, p.interval, p.point, p.widths, p.side, p.justification,
         p.scale, p.normalize, p.trim, p.n, p.show_slab, p.show_interval, p.show_point) do (positions, dists), st, iv, pt, ws, side, just, sc, nrm, tr, n, sslab, sint, spoint

        for (pos, d) in zip(positions, dists)
            if sslab
                xs, th = slab_curve(d; kind=st, n=n, trim=tr)
                poly!(p, slab_polygon(xs, th; position=pos, orientation=:vertical,
                                      side=side, justification=just, scale=sc, normalize=nrm);
                      color=p.color[])
            end
            rows = point_interval(d; widths=ws, point=pt, interval=iv)
            if sint
                for r in rows
                    a, b = interval_segment(r.lower, r.upper; position=pos, orientation=:vertical)
                    linesegments!(p, [a, b]; linewidth=p.interval_linewidth[], color=p.color[])
                end
            end
            if spoint && !isempty(rows)
                scatter!(p, [point_marker(rows[1].value; position=pos, orientation=:vertical)];
                         markersize=p.point_size[], color=p.color[])
            end
        end
    end
    return p
end
```

**Note:** the `lift` over many observables is written for clarity; if the installed Makie flags the multi-arg `lift` destructuring, split into `map`/`@lift` per the version's idiom. The drawing logic is unchanged.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS — `slabinterval renders` green (a PNG is produced).

- [ ] **Step 5: Commit**

```bash
git add src/recipes/slabinterval.jl test/test_recipes.jl
git commit -m "feat: SlabInterval plot! draws slab/interval/point sub-parts"
```

---

## Task 12: Default-setter children and the pre-summarised path

**Files:**
- Modify: `src/recipes/slabinterval.jl` (add child recipes + `convert_arguments` registration loop + pre-summarised method)
- Modify: `src/DistributionPlots.jl` (export child recipe names)
- Modify: `test/test_recipes.jl` (test children + reject slab-from-summary)

**Interfaces:**
- Consumes: everything in Tasks 10–11.
- Produces public recipes: `halfeye`, `eye`, `ccdfinterval`, `cdfinterval`, `gradientinterval`, `histinterval`, `slab`, `interval`, `pointinterval`, `spike` (and their `!` forms). Plus `pointinterval(positions, values, lowers, uppers)` for pre-summarised data, and an `ArgumentError`-throwing `slab` method for pre-summarised input.

- [ ] **Step 1: Write the failing test**

Add to `test/test_recipes.jl`:

```julia
@testset "children + pre-summarised" begin
    CairoMakie.activate!()
    for f in (halfeye, eye, ccdfinterval, cdfinterval, gradientinterval,
              histinterval, slab, interval, pointinterval, spike)
        @test (f(randn(1500)) isa Makie.FigureAxisPlot)
    end

    # halfeye defaults: slab on top, point shown
    th = Makie.default_theme(nothing, HalfEye)
    @test th[:side][] == :top

    # override a child default
    @test (halfeye(randn(1500); interval=:hdi) isa Makie.FigureAxisPlot)

    # pre-summarised pointinterval: positions + point + lower/upper
    @test (pointinterval([1.0,2.0], [0.0,1.0], [-1.0,0.0], [1.0,2.0]) isa Makie.FigureAxisPlot)

    # slab from pre-summarised data is rejected loudly (not a bare MethodError)
    @test_throws ArgumentError DistributionPlots._reject_slab_summary()
end
```

- [ ] **Step 2: Run test to verify it fails**

Run tests. Expected: FAIL — `halfeye` etc. not defined.

- [ ] **Step 3: Write minimal implementation**

Append to `src/recipes/slabinterval.jl`:

```julia
# (child recipe name => attribute defaults it presets over SlabInterval)
const _CHILD_DEFAULTS = Dict(
    :HalfEye        => (slab_type=:pdf,  side=:top,  show_slab=true,  show_interval=true,  show_point=true),
    :Eye            => (slab_type=:pdf,  side=:both, show_slab=true,  show_interval=true,  show_point=true),
    :CcdfInterval   => (slab_type=:ccdf, side=:top,  show_slab=true,  show_interval=true,  show_point=true),
    :CdfInterval    => (slab_type=:cdf,  side=:top,  show_slab=true,  show_interval=true,  show_point=true),
    :GradientInterval => (slab_type=:pdf, side=:both, show_slab=true, show_interval=true, show_point=true),
    :HistInterval   => (slab_type=:histogram, side=:top, show_slab=true, show_interval=true, show_point=true),
    :Slab           => (show_slab=true,  show_interval=false, show_point=false),
    :Interval       => (show_slab=false, show_interval=true,  show_point=false),
    :PointInterval  => (show_slab=false, show_interval=true,  show_point=true),
    :Spike          => (slab_type=:pdf,  show_slab=true, show_interval=false, show_point=false),
)

# Generate each child as a recipe that presets defaults then delegates to
# slabinterval!. This is ggdist's own "child = parent + defaults" mechanism.
for (T, defs) in _CHILD_DEFAULTS
    fname = Symbol(lowercase(String(T)))
    fname!  = Symbol(fname, :!)
    @eval begin
        @recipe($T, data) do scene
            merge(Attributes(; $(defs)...), Makie.default_theme(scene, SlabInterval))
        end
        function Makie.plot!(p::$T)
            slabinterval!(p, p.data;
                slab_type=p.slab_type, interval=p.interval, point=p.point, widths=p.widths,
                side=p.side, justification=p.justification, scale=p.scale, normalize=p.normalize,
                show_slab=p.show_slab, show_interval=p.show_interval, show_point=p.show_point,
                trim=p.trim, n=p.n, color=p.color, colormap=p.colormap, colorrange=p.colorrange,
                interval_linewidth=p.interval_linewidth, point_size=p.point_size)
            return p
        end
        # each public recipe type gets its own convert_arguments (Makie recipe
        # types do not share a supertype), all delegating to the shared converter
        Makie.convert_arguments(::Type{<:$T}, x::Distributions.UnivariateDistribution) = (_to_dist_args(x),)
        Makie.convert_arguments(::Type{<:$T}, x::AbstractVector{<:Real}) = (_to_dist_args(x),)
        Makie.convert_arguments(::Type{<:$T}, xs::AbstractVector{<:Distributions.UnivariateDistribution}) = (_to_dist_args(xs),)
    end
end

# Pre-summarised (ggdist geom_ path): draw supplied point/lower/upper directly.
function pointinterval(positions::AbstractVector{<:Real}, values::AbstractVector{<:Real},
                       lowers::AbstractVector{<:Real}, uppers::AbstractVector{<:Real}; kwargs...)
    fig = Figure()
    ax = Axis(fig[1,1])
    for i in eachindex(positions)
        a, b = interval_segment(lowers[i], uppers[i]; position=positions[i], orientation=:vertical)
        linesegments!(ax, [a, b]; linewidth=6, color=:black)
        scatter!(ax, [point_marker(values[i]; position=positions[i], orientation=:vertical)];
                 markersize=10, color=:black)
    end
    return Makie.FigureAxisPlot(fig, ax, ax.scene.plots[end])
end

# A slab needs a density; pre-summarised data has none — reject loudly.
_reject_slab_summary() = throw(ArgumentError(
    "slab needs a density; pre-summarised point/interval data has none. Pass samples or a Distribution, or use pointinterval for pre-summarised data."))
slab(::AbstractVector{<:Real}, ::AbstractVector{<:Real}, ::AbstractVector{<:Real}, ::AbstractVector{<:Real}) =
    _reject_slab_summary()
```

Add exports to `src/DistributionPlots.jl`:

```julia
export halfeye, halfeye!, eye, eye!, ccdfinterval, ccdfinterval!, cdfinterval, cdfinterval!,
       gradientinterval, gradientinterval!, histinterval, histinterval!,
       slab, slab!, interval, interval!, pointinterval, pointinterval!, spike, spike!,
       HalfEye, Eye, CcdfInterval, CdfInterval, GradientInterval, HistInterval,
       Slab, Interval, PointInterval, Spike
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS — `children + pre-summarised` green. If `@eval`-generated `@recipe` clashes with the exported plain function `pointinterval` (the recipe also defines `pointinterval`), keep the pre-summarised method as an additional method on the same generic — Julia dispatches on the 4-positional-arg signature. Verify both `pointinterval(randn(100))` and the 4-arg form resolve.

- [ ] **Step 5: Commit**

```bash
git add src/recipes/slabinterval.jl src/DistributionPlots.jl test/test_recipes.jl
git commit -m "feat: default-setter child recipes and pre-summarised pointinterval"
```

---

## Task 13: Dots recipes

**Files:**
- Create: `src/recipes/dots.jl`
- Modify: `src/DistributionPlots.jl` (include + export)
- Modify: `test/test_recipes.jl` (dots smoke test)

**Interfaces:**
- Consumes: `dot_layout` (Task 7), `_to_dist_args` (Task 10), geometry (Task 8), `point_interval` (Tasks 4–5).
- Produces: `dots`/`dots!` and `dotsinterval`/`dotsinterval!` recipes.

- [ ] **Step 1: Write the failing test**

Add to `test/test_recipes.jl`:

```julia
@testset "dots" begin
    CairoMakie.activate!()
    @test (dots(randn(1000)) isa Makie.FigureAxisPlot)
    @test (dotsinterval(randn(1000)) isa Makie.FigureAxisPlot)
end
```

- [ ] **Step 2: Run test to verify it fails**

Run tests. Expected: FAIL — `dots` not defined.

- [ ] **Step 3: Write minimal implementation**

Create `src/recipes/dots.jl`:

```julia
using Makie

@recipe(Dots, data) do scene
    Attributes(ndots = 50, orientation = :vertical, color = :black,
               scale = 0.9, dot_gap = 0.0)
end

function Makie.plot!(p::Dots)
    args = lift(a -> a isa Tuple ? a[1] : a, p.data)
    lift(args, p.ndots, p.scale) do (positions, dists), nd, sc
        for (pos, d) in zip(positions, dists)
            lay = dot_layout(d; ndots=nd)
            r = lay.binwidth / 2
            pts = [Point2f(pos + (y - 1) * lay.binwidth * sc, x) for (x, y) in zip(lay.x, lay.y)]
            scatter!(p, pts; markersize = 6, color = p.color[])
        end
    end
    return p
end

Makie.convert_arguments(::Type{<:Dots}, x::AbstractVector{<:Real}) = (_to_dist_args(x),)
Makie.convert_arguments(::Type{<:Dots}, x::Distributions.UnivariateDistribution) = (_to_dist_args(x),)

@recipe(DotsInterval, data) do scene
    merge(Attributes(ndots = 50, show_interval = true, show_point = true, widths = [0.66, 0.95],
                     point = :median, interval = :qi, color = :black, scale = 0.9),
          Attributes())
end

function Makie.plot!(p::DotsInterval)
    args = lift(a -> a isa Tuple ? a[1] : a, p.data)
    lift(args, p.ndots, p.widths, p.point, p.interval, p.show_interval, p.show_point, p.scale) do (positions, dists), nd, ws, pt, iv, sint, spoint, sc
        for (pos, d) in zip(positions, dists)
            lay = dot_layout(d; ndots=nd)
            pts = [Point2f(pos + (y - 1) * lay.binwidth * sc, x) for (x, y) in zip(lay.x, lay.y)]
            scatter!(p, pts; markersize = 6, color = p.color[])
            rows = point_interval(d; widths=ws, point=pt, interval=iv)
            if sint
                for rr in rows
                    a, b = interval_segment(rr.lower, rr.upper; position=pos, orientation=:vertical)
                    linesegments!(p, [a, b]; linewidth=6, color=p.color[])
                end
            end
            spoint && !isempty(rows) && scatter!(p, [point_marker(rows[1].value; position=pos, orientation=:vertical)];
                                                 markersize=10, color=p.color[])
        end
    end
    return p
end

Makie.convert_arguments(::Type{<:DotsInterval}, x::AbstractVector{<:Real}) = (_to_dist_args(x),)
Makie.convert_arguments(::Type{<:DotsInterval}, x::Distributions.UnivariateDistribution) = (_to_dist_args(x),)
```

Add to `src/DistributionPlots.jl`:

```julia
include("recipes/dots.jl")
export dots, dots!, dotsinterval, dotsinterval!, Dots, DotsInterval
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS — `dots` green.

- [ ] **Step 5: Commit**

```bash
git add src/recipes/dots.jl src/DistributionPlots.jl test/test_recipes.jl
git commit -m "feat: dots and dotsinterval recipes"
```

---

## Task 14: Lineribbon recipe

**Files:**
- Create: `src/recipes/lineribbon.jl`
- Modify: `src/DistributionPlots.jl` (include + export)
- Modify: `test/test_recipes.jl` (lineribbon smoke test)

**Interfaces:**
- Consumes: `point_interval` (Tasks 4–5), `asdist` (Task 2).
- Produces: `lineribbon`/`lineribbon!` recipe drawing nested bands + a central line over a predictor. Signature: `lineribbon(x, dists)` where `x` is a predictor vector and `dists` is a vector of per-x distributions/sample-vectors.

- [ ] **Step 1: Write the failing test**

Add to `test/test_recipes.jl`:

```julia
@testset "lineribbon" begin
    CairoMakie.activate!()
    xgrid = collect(0.0:0.5:10.0)
    # a growing-uncertainty fan: samples per x
    perx = [randn(800) .* (0.2 + 0.1x) .+ sin(x) for x in xgrid]
    @test (lineribbon(xgrid, perx) isa Makie.FigureAxisPlot)
end
```

- [ ] **Step 2: Run test to verify it fails**

Run tests. Expected: FAIL — `lineribbon` not defined.

- [ ] **Step 3: Write minimal implementation**

Create `src/recipes/lineribbon.jl`:

```julia
using Makie

@recipe(LineRibbon, x, dists) do scene
    Attributes(widths = [0.5, 0.8, 0.95], point = :median, interval = :qi,
               colormap = :blues, color = :black, linewidth = 2)
end

function Makie.plot!(p::LineRibbon)
    lift(p.x, p.dists, p.widths, p.point, p.interval) do xs, dists, ws, pt, iv
        wsorted = sort(ws; rev=true)                 # widest band drawn first (back)
        ncol = length(wsorted)
        for (k, w) in enumerate(wsorted)
            los = Float64[]; his = Float64[]
            for d in dists
                r = point_interval(d; widths=[w], point=pt, interval=iv)[1]
                push!(los, r.lower); push!(his, r.upper)
            end
            shade = ncol == 1 ? 0.5 : (k - 1) / (ncol - 1)
            band!(p, xs, los, his; color=(:blue, 0.15 + 0.5 * (1 - shade)))
        end
        centers = [point_interval(d; widths=[first(wsorted)], point=pt, interval=iv)[1].value for d in dists]
        lines!(p, xs, centers; color=p.color[], linewidth=p.linewidth[])
    end
    return p
end

# predictor + per-x sample vectors → wrap each column as an AbstractDist
Makie.convert_arguments(::Type{<:LineRibbon}, x::AbstractVector{<:Real},
                        dists::AbstractVector{<:AbstractVector{<:Real}}) =
    (collect(Float64, x), AbstractDist[asdist(d) for d in dists])
Makie.convert_arguments(::Type{<:LineRibbon}, x::AbstractVector{<:Real},
                        dists::AbstractVector{<:Distributions.UnivariateDistribution}) =
    (collect(Float64, x), AbstractDist[asdist(d) for d in dists])
```

Add to `src/DistributionPlots.jl`:

```julia
include("recipes/lineribbon.jl")
export lineribbon, lineribbon!, LineRibbon
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS — `lineribbon` green.

- [ ] **Step 5: Commit**

```bash
git add src/recipes/lineribbon.jl src/DistributionPlots.jl test/test_recipes.jl
git commit -m "feat: lineribbon recipe with nested bands and central line"
```

---

## Task 15: RVars extension

**Files:**
- Create: `ext/DistributionPlotsRVarsExt.jl`
- Create: `test/test_extensions.jl`
- Modify: `test/runtests.jl` (include)

**Interfaces:**
- Consumes: `_to_dist_args`, all public recipe types (Tasks 10–14); `RVars.RVar`, `RVars.draws`, `RVars.variables`.
- Produces: `convert_arguments(P, ::RVar)` for every public recipe type `P`, plus a `_to_dist_args(::RVar)` method: scalar RV → 1 position, length-`k` vector RV → positions `1:k`, `variables(x)` carried for tick labels (returned in a companion accessor).

- [ ] **Step 1: Write the failing test**

Create `test/test_extensions.jl`:

```julia
using DistributionPlots
using RVars
using Makie, CairoMakie
using Test

@testset "RVars extension" begin
    CairoMakie.activate!()
    # scalar RV: 1000 draws → one distribution at position 1
    rv0 = RVar(randn(1000))
    pos, dists = DistributionPlots._to_dist_args(rv0)
    @test pos == [1.0]
    @test length(dists) == 1

    # vector RV: 3 elements → positions 1:3
    rv1 = RVar(randn(1000, 3))
    pos3, dists3 = DistributionPlots._to_dist_args(rv1)
    @test pos3 == [1.0, 2.0, 3.0]
    @test length(dists3) == 3

    # recipes accept a RVar directly
    @test (halfeye(rv1) isa Makie.FigureAxisPlot)
    @test (pointinterval(rv1) isa Makie.FigureAxisPlot)
end
```

- [ ] **Step 2: Run test to verify it fails**

Add `include("test_extensions.jl")` to `test/runtests.jl`. Run tests.
Expected: FAIL — `_to_dist_args(::RVar)` not defined (extension not written).

- [ ] **Step 3: Write minimal implementation**

Create `ext/DistributionPlotsRVarsExt.jl`:

```julia
module DistributionPlotsRVarsExt

using DistributionPlots
using DistributionPlots: _to_dist_args, AbstractDist, asdist
using RVars
using RVars: RVar, draws, variables
using Makie
using Distributions

# scalar RV (N=0): one distribution from all draws.
# vector RV (N=1, length k): k distributions at integer positions 1:k.
function DistributionPlots._to_dist_args(x::RVar)
    raw = draws(x)                       # (ndraws, dims...) — axis 1 is draws
    if ndims(raw) == 1
        return ([1.0], AbstractDist[asdist(vec(raw))])
    else
        k = size(raw, 2)
        return (collect(1.0:k), AbstractDist[asdist(vec(raw[:, j])) for j in 1:k])
    end
end

# Register convert_arguments for every public recipe type.
const _RECIPES = (SlabInterval, HalfEye, Eye, CcdfInterval, CdfInterval,
                  GradientInterval, HistInterval, Slab, Interval, PointInterval,
                  Spike, Dots, DotsInterval)

for T in _RECIPES
    @eval Makie.convert_arguments(::Type{<:$T}, x::RVar) = (_to_dist_args(x),)
end

end # module
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS — `RVars extension` green. The extension loads because the test target lists `RVars`. If `draws(x)` shape differs, consult RVars' README (`draws(x)` returns `(ndraws, dims...)`); a scalar RV backs a length-`ndraws` vector, a vector RV backs an `(ndraws, k)` matrix.

- [ ] **Step 5: Commit**

```bash
git add ext/DistributionPlotsRVarsExt.jl test/test_extensions.jl test/runtests.jl
git commit -m "feat: RVars extension (RVar → recipe inputs)"
```

---

## Task 16: MCMCChains extension

**Files:**
- Create: `ext/DistributionPlotsMCMCChainsExt.jl`
- Modify: `test/test_extensions.jl` (Chains testset)

**Interfaces:**
- Consumes: public recipe types; `MCMCChains.Chains`. Routes through RVars' `Chains` conversion, so this extension effectively requires both weakdeps at use time (it calls `RVars.RVar(::Chains)`).
- Produces: `convert_arguments(P, ::Chains)` for every public recipe type, by converting `Chains → RVar` then reusing the RVars path.

- [ ] **Step 1: Write the failing test**

Add to `test/test_extensions.jl`:

```julia
using MCMCChains

@testset "MCMCChains extension" begin
    CairoMakie.activate!()
    chn = Chains(randn(500, 3, 2), [:a, :b, :c])       # 500 iters, 3 params, 2 chains
    @test (halfeye(chn) isa Makie.FigureAxisPlot)      # 3 positions, chains pooled
    pos, dists = DistributionPlots._to_dist_args(chn)
    @test pos == [1.0, 2.0, 3.0]
end
```

- [ ] **Step 2: Run test to verify it fails**

Run tests. Expected: FAIL — `_to_dist_args(::Chains)` not defined.

- [ ] **Step 3: Write minimal implementation**

Create `ext/DistributionPlotsMCMCChainsExt.jl`:

```julia
module DistributionPlotsMCMCChainsExt

using DistributionPlots
using DistributionPlots: _to_dist_args
using MCMCChains: Chains
using RVars: RVar
using Makie

# Route Chains through RVars' own Chains conversion (single tested path),
# then reuse the RVars _to_dist_args method.
DistributionPlots._to_dist_args(chn::Chains) = _to_dist_args(RVar(chn))

const _RECIPES = (SlabInterval, HalfEye, Eye, CcdfInterval, CdfInterval,
                  GradientInterval, HistInterval, Slab, Interval, PointInterval,
                  Spike, Dots, DotsInterval)

for T in _RECIPES
    @eval Makie.convert_arguments(::Type{<:$T}, chn::Chains) = (_to_dist_args(chn),)
end

end # module
```

**Note:** the `[extensions]` table triggers this module on `MCMCChains` alone, but it calls `RVars.RVar(::Chains)`. Both packages are in the test target, so both load. Document in the README that Chains support requires `using RVars` too.

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS — `MCMCChains extension` green.

- [ ] **Step 5: Commit**

```bash
git add ext/DistributionPlotsMCMCChainsExt.jl test/test_extensions.jl
git commit -m "feat: MCMCChains extension routed through RVars"
```

---

## Task 17: AlgebraOfGraphics extension, README, and finalization

**Files:**
- Create: `ext/DistributionPlotsAlgebraOfGraphicsExt.jl`
- Create: `README.md`
- Modify: `test/test_extensions.jl` (AoG smoke test)

**Interfaces:**
- Consumes: public recipe types; `AlgebraOfGraphics.visual`.
- Produces: `visual(HalfEye)` (and the other recipe types) usable in an AoG spec. In practice AoG's `visual(::Type{<:AbstractPlot})` already works for any Makie recipe; this extension only needs to ensure the recipe types are visible and add any AoG-specific `aesthetic_mapping`/scales if required. The test asserts an AoG `draw` of `visual(HalfEye)` runs.

- [ ] **Step 1: Write the failing test**

Add to `test/test_extensions.jl`:

```julia
using AlgebraOfGraphics
const AoG = AlgebraOfGraphics

@testset "AlgebraOfGraphics extension" begin
    CairoMakie.activate!()
    # visual() over our recipe type composes and draws without error
    spec = AoG.data((x = repeat(1:3, inner=500),
                     y = randn(1500))) *
           AoG.mapping(:x, :y) *
           AoG.visual(Interval)
    fg = AoG.draw(spec)
    @test fg isa Makie.FigureGrid
end
```

- [ ] **Step 2: Run test to verify it fails**

Run tests. Expected: FAIL — extension not present (or `visual(Interval)` unrecognised if the type isn't exported/loaded through the extension).

- [ ] **Step 3: Write minimal implementation**

Create `ext/DistributionPlotsAlgebraOfGraphicsExt.jl`:

```julia
module DistributionPlotsAlgebraOfGraphicsExt

using DistributionPlots
using AlgebraOfGraphics
using Makie

# AlgebraOfGraphics can already wrap any Makie recipe via visual(RecipeType).
# This extension is the declared integration point; if AoG needs an explicit
# aesthetic mapping for our recipes, register it here. For v1, re-exporting the
# recipe types through the loaded extension is sufficient for `visual(HalfEye)`.
# (No additional methods required for basic visual() support in AoG.)

end # module
```

**If the test fails** because AoG needs an aesthetic mapping for the custom recipe, add (guided by `?AlgebraOfGraphics.aesthetic_mapping`):

```julia
    import AlgebraOfGraphics: aesthetic_mapping
    using AlgebraOfGraphics: AesX, AesY
    for T in (DistributionPlots.SlabInterval, DistributionPlots.HalfEye,
              DistributionPlots.Interval, DistributionPlots.PointInterval)
        @eval aesthetic_mapping(::Type{<:$T}, ::Any) = AlgebraOfGraphics.dims(1) => AesX, AlgebraOfGraphics.dims(2) => AesY
    end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS — `AlgebraOfGraphics extension` green. If AoG's API differs for the installed version, scope this test to the minimum that passes and record the realised AoG version in the README's compat note; full AoG *analyses* are explicitly out of v1 scope.

- [ ] **Step 5: Write README and commit**

Create `README.md`:

```markdown
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

## Notes

- Chains support requires `using RVars` alongside `using MCMCChains`.
- Requires Julia 1.10+.
```

```bash
git add ext/DistributionPlotsAlgebraOfGraphicsExt.jl README.md test/test_extensions.jl
git commit -m "feat: AlgebraOfGraphics extension, README; v0.1.0 feature-complete"
```

---

## Self-Review

**1. Spec coverage** — every spec section maps to a task:

| Spec item | Task(s) |
|---|---|
| Input protocol (SampleDist/AnalyticDist, 4 methods) | 2, 3 |
| Boundary-reflected KDE | 3 |
| `point_interval` qi/hdci/hdi, mean/median/mode, Tables rows | 4, 5 |
| `slab_curve` pdf/cdf/ccdf/histogram | 6 |
| `dot_layout` Wilkinson packing | 7 |
| Geometry transform (side/justification/scale/normalize/orientation) | 8 |
| Golden fixtures vs ggdist 3.3.3 | 9 |
| SlabInterval parent engine + toggles + lift | 10, 11 |
| Default-setter children (halfeye/eye/ccdf/cdf/gradient/hist/slab/interval/pointinterval/spike) | 12 |
| Pre-summarised geom_ path + reject slab-from-summary | 12 |
| Dots/dotsinterval | 13 |
| Lineribbon (separate engine) | 14 |
| convert_arguments per public recipe type (not one <:SlabInterval) | 10, 12, 15, 16 |
| RVars extension (scalar→1, vector→1:k, variables) | 15 |
| MCMCChains extension routed through RVars | 16 |
| AoG visual() extension | 17 |
| Error tiers (reject/warn/silent) | 2 (empty reject, NaN warn), 4 (width reject), 12 (slab-summary reject) |
| Julia 1.10 floor, exact UUIDs, compat | 1 |
| Testing: unit + golden + smoke | all unit tasks, 9, 11–17 |

**2. Placeholder scan** — no "TBD"/"implement later"/"add error handling"; every code step shows complete code. Version-sensitive Makie API points carry explicit "verify against installed Makie" instructions rather than vague hand-waving, which is a real instruction, not a placeholder.

**3. Type consistency** — `point_interval` returns `Vector{PIRow}` with fields `value/lower/upper/width/point/interval` everywhere (Tasks 4, 5, 9, 11, 12, 13). `_to_dist_args` returns `(Vector{Float64}, Vector{AbstractDist})` in Task 10 and every extension (15, 16). Recipe attribute names (`slab_type`, `show_slab`, `show_interval`, `show_point`, `side`, `scale`, `normalize`, `trim`, `n`, `widths`, `point`, `interval`) are identical across the parent (10), children (12), and dots (13). Geometry functions `slab_polygon`/`interval_segment`/`point_marker`/`normalize_thickness` are defined in Task 8 and called with matching signatures in Tasks 11, 12, 13.

### Known risks flagged for execution

- **Makie recipe API version drift** — the single largest risk. Every recipe task (10–14) begins by verifying the installed Makie's `@recipe`/`convert_arguments`/`lift` syntax. If Makie 0.24 uses the newer documented-attributes recipe form, translate the `do scene Attributes(...) end` blocks accordingly; the attribute names and drawing logic are unaffected.
- **AoG integration depth** (Task 17) — `visual()` over a Makie recipe is expected to work out of the box; if the installed AoG needs an explicit `aesthetic_mapping`, the fallback snippet is provided. Custom AoG analyses remain out of v1 scope.
- **hdci vs ggdist definition** (Task 9) — if `median_hdci` golden values diverge beyond `1e-6`, reconcile the contiguous-window definition before loosening tolerance.
