# PlotDist.jl — Design Specification

**Date:** 2026-07-24
**Status:** Approved design, pre-implementation
**Author:** Karim Naguib (with Claude Code)

## Purpose

`PlotDist.jl` provides Makie plots and recipes that mimic R's
[`ggdist`](https://mjskay.github.io/ggdist/) — slab, interval, pointinterval,
lineribbon, and the dots/dotsinterval family — for visualizing distributions
and uncertainty. It plots both **sample-based** random variables (primarily
[`RandomDraws.jl`](https://github.com/karimn/RandomDraws.jl)'s `RandomDraw`, the
Julia analog of `posterior::rvar`) and **analytic** distributions
(`Distributions.jl`), through one interface, the way ggdist accepts both `rvar`
and `distributional` objects.

The package is faithful to ggdist's **class hierarchy** and to **Makie's recipe
conventions** simultaneously; where the two frameworks differ, the deviations
are deliberate and documented below.

## Scope

### In scope (v1)

- **Slab/interval family**, mirroring ggdist's hierarchy: a `slabinterval`
  parent engine with `slab`, `interval`, `pointinterval`, `dotsinterval`,
  `dots`, and `spike` as specializations, and `halfeye`, `eye`, `ccdfinterval`,
  `cdfinterval`, `gradientinterval`, `histinterval` as thin default-setters over
  the parent.
- **`lineribbon`** — a separate engine (a sibling of the slabinterval family in
  ggdist), drawing nested probability bands over a predictor with a central
  line.
- **Point/interval summaries** — `qi` (equal-tailed quantile interval), `hdci`
  (highest-density continuous interval), `hdi` (highest-density interval, may be
  disjoint); point summaries `mean`, `median`, `mode`.
- **Four input types** through one interface: `RandomDraw`, raw
  `AbstractVector`/`AbstractMatrix` of samples, `Distributions.Distribution`,
  and `MCMCChains.Chains`. Plus **pre-summarised** point/lower/upper data
  (ggdist's `geom_` path).
- **A public `point_interval` function** returning a `Tables.jl`-compatible row
  set — a tidybayes-style summary usable without plotting.
- **AlgebraOfGraphics** compatibility via `visual()` (extension).

### Out of scope (v1)

- Custom AlgebraOfGraphics *analyses* (so that `mapping(...) *
  pointinterval_analysis()` composes) — stretch goal, not v1.
- Pixel-level visual regression testing.
- Full discrete-distribution slab support beyond straightforward pmf cases;
  unsupported discrete cases error clearly rather than draw something wrong.
- Per-chain plotting (chains are pooled in v1).

## Architecture

### The reusable foundation is two stat functions, not a shared geom

```
                    ┌─ point_interval(d, widths; point, interval) ─┐
  input ─► d ┤ (normalised to SampleDist or AnalyticDist)          ├─► recipes
                    └─ slab_curve(d; kind, n) ─────────────────────┘
                       dot_layout(d; n)

  point_interval  consumed by  interval, pointinterval, lineribbon
  slab_curve      consumed by  slab, spike
  dot_layout      consumed by  dots, dotsinterval
```

Makie has **no native stat/geom separation** (unlike ggplot2 — see "Makie
conventions" below), so this separation is a **convention we enforce** through
code layout and testing, not a framework feature. The pure stat layer is what
makes both the R golden-fixture tests and Makie's reactive `lift` wrapping
possible.

### Three internal layers (none are recipes)

| Layer | Files | Responsibility |
|---|---|---|
| Input protocol | `interface.jl` | The `SampleDist` / `AnalyticDist` types and their four-method interface. Exposed to Makie via `convert_arguments`. |
| Statistics | `intervals.jl`, `slabs.jl`, `density.jl`, `dotlayout.jl` | Pure functions → plain numbers. Point/interval summaries, slab curves, KDE (incl. boundary reflection), Wilkinson dot packing. |
| Geometry | `geometry.jl` | `(thickness, position, side, justification, scale, orientation, normalize)` → drawing coordinates. Applied identically to every sub-part. |

### The input protocol (four methods, two types)

Every accepted input normalises — **at the `convert_arguments` boundary** — to
one of two internal types, each answering the same four questions:

| Method | `SampleDist` (RandomDraw / Vector / Chains) | `AnalyticDist` (Distributions.jl) |
|---|---|---|
| `support(d)` | `(minimum, maximum)` of draws | `Distributions.support`, trimmed for infinite tails to a quantile cutoff |
| `quantile_at(d, p)` | `StatsBase.quantile` (type-7, matches R's default) | `Distributions.quantile` (exact) |
| `cdf_at(d, x)` | empirical CDF | `Distributions.cdf` (exact) |
| `density_at(d, xs)` | KDE on `xs` (boundary-reflected when support is bounded) | `Distributions.pdf` (exact) |

This mirrors modern ggdist, which unified its `x` (samples) and `xdist`
(distribution) paths into one stat that dispatches on which is supplied — hence
our two types feed *one* engine, not two.

## Class hierarchy (faithful to ggdist 3.3.3)

Verified against ggdist source. Geom tree:

```
AbstractGeom
├── GeomSlabinterval          ← the real primitive
│   ├── GeomSlab
│   │   └── GeomSpike
│   ├── GeomInterval
│   ├── GeomPointinterval
│   └── GeomDotsinterval
│       └── GeomDots
└── GeomLineribbon            ← sibling, outside the family
```

**Key finding:** ggdist's `stat_halfeye`, `stat_eye`, `stat_ccdfinterval`,
`stat_cdfinterval`, `stat_gradientinterval`, `stat_histinterval` are **not
distinct classes** — they are `stat_slabinterval` with preset defaults,
manufactured by a factory. Likewise the 14 `stat_dist_*` exports are the
analytic-input path of the same stats. So ggdist's ~28 `stat_` exports collapse
to a handful of real classes plus defaults. We reproduce exactly that:
convenience variants are default-setters, never new drawing logic. Building them
as independent primitives would *diverge* from ggdist.

### How the hierarchy maps to Julia

| ggdist | Julia mirror | Faithful? |
|---|---|---|
| Stat ggproto tree | **Abstract type tree** (`AbstractSlabinterval`, `AbstractDotsinterval`, …) for computation dispatch | ✅ exact |
| `stat_halfeye` = defaults over `stat_slabinterval` | `halfeye!` calls `slabinterval!` with preset attributes | ✅ exact (same mechanism) |
| Geom inheritance (`GeomSlab <: GeomSlabinterval`) | Makie `@recipe` generates concrete `Plot{f}` types that can't subtype each other | ⚠️ approximated via parent-engine + default-children |

The one gap — Makie recipe *types* can't inherit — has no user-visible
consequence. The hierarchy is faithfully reproduced everywhere it carries
meaning: real `<:` inheritance in the stat layer, and ggdist's own
parent+defaults mechanism in the recipe layer.

## Recipe layer: "structure like A, use like C"

One `slabinterval` engine is the parent; the pieces are its specializations;
each piece is still independently callable and stackable — exactly as ggdist's
`geom_slab` and `geom_interval` are both children of `slabinterval` *and*
usable on their own.

```
slabinterval!  ── parent engine (draws any enabled sub-part)
  ├── slab!            (slab sub-part only)
  │   └── spike!       (slab collapsed to one x)
  ├── interval!        (interval sub-part only)
  ├── pointinterval!   (interval + point)
  ├── dotsinterval!    (dots + interval)
  │   └── dots!        (dots sub-part only)
  └── default-setters over slabinterval!:
        halfeye! eye! ccdfinterval! cdfinterval! gradientinterval! histinterval!

lineribbon!    ── separate engine (nested bands + central line over a predictor)
```

### The parent engine

`slabinterval!` does three things per call:

1. Normalise input → `SampleDist`/`AnalyticDist` (via `convert_arguments`).
2. Call the stat functions its attributes request, **wrapped in `lift`** so the
   plot is reactive.
3. Draw whichever sub-parts are enabled, all sharing one positioning transform.

- **Sub-part toggles** mirror ggdist's `show_slab`/`show_point`/`show_interval`.
  A child presets them: `slab!` = slab on, others off; `pointinterval!` =
  interval + point on; `halfeye!` = slab + interval + point with `side=:top`.
- **Shared positioning contract** — `orientation ∈ {:vertical,:horizontal}`,
  `side ∈ {:top,:bottom,:both}`, `justification`, `scale`, `normalize` — applied
  by one `geometry.jl` transform to *every* sub-part, so slab and interval in
  the same call align by construction.
- **Children are default-setters, not new logic.** For example
  `halfeye!(ax, x, d)` ≡ `slabinterval!(ax, x, d; slab_type=:pdf, side=:top,
  show_point=true, interval=:qi)`. New variant = new default set (~3 lines).
- Users can override a child's defaults (`halfeye!(...; interval=:hdi)`) because
  `slab_type`/`interval`/`point` are **attributes**, not hardcoded.
- **Composite convenience recipes remain full Makie recipes** (single plot
  object → legend entry, AoG `visual()`-drivable), not plain functions.

## Makie conventions

The design matches Makie idiom, with three refinements that make it properly
idiomatic:

**Already conventional:**
- Paired lowercase `foo`/`foo!` recipe functions via `@recipe` (cf.
  `violin`/`violin!`).
- `plot!(ax, positions, data)` signature with `(position, values)` ordering
  (cf. `boxplot(xs, ys)`).
- Parent composing children via `plot!` calling other `plot!` — Makie's "full
  recipe" model (cf. `rainclouds` calling `violin!`/`boxplot!`/`scatter!`).

**Refinements folded in:**
1. **Input normalisation goes through `convert_arguments`** — Makie's designated
   hook for accepting custom input types — not a bespoke pre-step. Core defines
   it for `Distribution` and `AbstractVector`; extensions add `RandomDraw` and
   `Chains`.
2. **Stat computation is wrapped in `lift`/Observables** so the recipe is
   reactive (recomputes when inputs or attributes change). The pure stat
   functions are trivially `lift`-wrappable *because* they're pure.
3. **Interval-level → colour via `colormap`/`colorrange` attributes**
   (themeable, cycle-compatible), not a private palette.

**Attribute naming:** ggdist-style underscored names (`slab_type`,
`show_interval`, `point_interval`) are kept. Makie's *primitive* recipes avoid
underscores, but its *higher-level statistical* recipes (e.g. `rainclouds` with
`plot_boxplots`) use them — our tier matches that convention, so ggdist
familiarity and Makie idiom align with no trade-off. (Exact `rainclouds`
spellings to be confirmed against installed Makie during implementation.)

## Statistics layer detail

### `point_interval(d, widths; point, interval)`

The workhorse feeding `interval`, `pointinterval`, **and** `lineribbon`.

- `point ∈ {mean, median, mode}`. `mode` requires the KDE (argmax of density) —
  which is why density estimation is load-bearing for *reported numbers*, not
  just slab shape.
- `interval ∈ {qi, hdci, hdi}`:
  - `qi` — equal-tailed quantiles. Default. One interval.
  - `hdci` — narrowest contiguous window over sorted draws holding the target
    mass. One interval, no density needed.
  - `hdi` — highest-density region from the KDE level set. **May return several
    disjoint intervals** for multimodal distributions.
- **Returns a `Tables.jl`-compatible row set**: one row per (position × width ×
  interval-piece), columns `value, .lower, .upper, .width, .point, .interval`.
  The set-of-intervals shape is baked in here so `hdi` disjoint intervals flow
  to both the interval bar and the ribbon band without per-recipe special
  casing. Callable standalone as a tidybayes-style summary.

### `slab_curve(d; kind, n)`

Feeds `slab`/`spike` only. `kind ∈ {pdf, cdf, ccdf, histogram}`; returns
`(xs, thickness)`. `pdf` from KDE (samples) or exact `pdf` (analytic);
`histogram` bypasses KDE.

### `dot_layout(d; n, quantiles)`

Wilkinson dot-packing → per-dot `(x, y, radius)`. The one non-trivial algorithm;
ported directly and pinned against R fixtures.

## Error handling

### Three-tier philosophy

The sorting rule is *whether the result would mislead*.

| Tier | Rule | Mechanism |
|---|---|---|
| **Reject loudly** | Structurally invalid | `throw(ArgumentError)` at the `convert_arguments` boundary, offending value in the message |
| **Degrade with a warning** | Drawable but compromised, or a fallback substituted | `@warn` once, then proceed — never silent |
| **Handle silently** | Valid, just degenerate | Draw the sensible thing |

**Reject:** empty samples; `lower > upper` in pre-summarised input;
`width ∉ (0,1)`; a slab requested from pre-summarised input (no density exists —
enforced by dispatch, so it surfaces as a `MethodError`-class rejection).

**Degrade + warn:** `NaN`/`missing` in samples → drop and `@warn` the count
(matches ggplot `na.rm`; silently dropping would alter the posterior); `hdi` on
a sample too small for a stable density → compute but warn; any unsupported
analytic case that falls back → warn what was substituted.

**Handle silently:** constant / zero-variance sample → draw a spike at the value
(detect variance first so `kde` doesn't throw); a single HDI interval where
multimodality was possible → just one interval.

### Two correctness landmines (live in the stat layer, not drawing)

- **Infinite analytic support** (`Normal`, `Cauchy`: `minimum(d) == -Inf`) →
  trim the slab range to `quantile(d, trim)` / `quantile(d, 1-trim)`, with
  `trim` a themeable attribute defaulting to ggdist's `.001`.
- **Bounded-support samples** (e.g. proportions in `[0,1]`) → naive KDE leaks
  past the boundary and drags `hdi` endpoints inward. Matching ggdist's
  `density="bounded"` requires **boundary reflection**, which `KernelDensity.jl`
  does not do natively — we implement it in `slabs.jl`, pinned by a Beta-
  distribution golden fixture.

## Dependencies

**Core (hard deps):** `Makie`, `Statistics`, `StatsBase`, `KernelDensity`,
`Distributions`, `Tables`.

- `Distributions` is a direct dep despite being weakdep-shaped, because
  `KernelDensity` already pulls it in transitively — declaring it explicitly
  costs nothing and avoids an extension that would load unconditionally anyway.
- `Tables` (not `DataFrames`) is the output contract for `point_interval` — an
  interface package with a near-zero transitive tree, letting results flow into
  DataFrames/AoG/CSV without a heavy dependency.

**Extensions (weakdeps), each a thin, independent adapter:**

| Extension | Weakdep | Provides | If not loaded |
|---|---|---|---|
| `PlotDistRandomDrawsExt` | `RandomDraws` | `convert_arguments(::Type{<:SlabInterval}, ::RandomDraw)`; scalar RV → 1 position, length-`k` vector RV → `k` positions; `variables(x)` → tick labels | core still works on `Vector`/`Distribution` |
| `PlotDistMCMCChainsExt` | `MCMCChains` | `convert_arguments(..., ::Chains)`, routed through RandomDraws' existing `Chains` conversion (not reimplemented) | `Chains` not plottable |
| `PlotDistAlgebraOfGraphicsExt` | `AlgebraOfGraphics` | `visual(HalfEye)` compatibility | recipes still work directly |

Extensions are kept **separate** (not one combined `PlotDistStatsExt`) so that
loading `MCMCChains` without `RandomDraws` — or vice versa — can't break the
other. Chains support routing through RandomDraws means a Chains-axis-order bug
(the class of bug RandomDraws' own CLAUDE.md documents) can only exist in one
place.

## File layout (planned)

```
PlotDist.jl/
├── Project.toml
├── README.md
├── src/
│   ├── PlotDist.jl          # module, includes, exports
│   ├── interface.jl         # SampleDist / AnalyticDist, convert_arguments
│   ├── intervals.jl         # point_interval (qi/hdci/hdi), point summaries
│   ├── density.jl           # KDE incl. boundary reflection
│   ├── slabs.jl             # slab_curve
│   ├── dotlayout.jl         # Wilkinson dot packing
│   ├── geometry.jl          # positioning transform
│   └── recipes/
│       ├── slabinterval.jl  # parent engine + default-setter children
│       ├── dots.jl          # dotsinterval / dots / spike
│       └── lineribbon.jl    # separate engine
├── ext/
│   ├── PlotDistRandomDrawsExt/
│   ├── PlotDistMCMCChainsExt/
│   └── PlotDistAlgebraOfGraphicsExt/
└── test/
    ├── runtests.jl
    ├── fixtures/            # committed R-exported golden data
    ├── gen_fixtures.R       # re-runnable, pins ggdist/posterior versions
    └── ...
```

## Testing strategy

Three layers, matching the choices made during design:

1. **Unit tests on the pure stat layer** — `point_interval` (qi/hdci/hdi),
   `slab_curve`, `dot_layout`, the geometry transform, and each
   `convert_arguments` method. Numbers in, numbers out, no Makie.
2. **Golden fixtures from R ggdist 3.3.3 / posterior 1.6.1** — a committed,
   re-runnable R script (`gen_fixtures.R`) exports reference outputs to
   `test/fixtures/`; Julia asserts within per-case tolerance.
3. **Smoke tests** — every recipe constructs and draws without error on a
   headless CairoMakie backend.

**Same-input guarantee:** sample-based fixtures export *both* the raw draws and
ggdist's output. Julia loads ggdist's *exact* input draws — never regenerated
"similar" data — so a mismatch means a real algorithm bug, not a cross-language
RNG difference.

**Fixture cases chosen to break things:**

| Distribution | Catches |
|---|---|
| `Normal` | symmetric baseline — `qi == hdi`, sanity floor |
| `Beta(2,8)` | bounded + asymmetric — boundary-reflection KDE; wrong bandwidth drifts `hdi` |
| bimodal mixture | disjoint `hdi` — proves the set-of-intervals return type |
| `StudentT` heavy tails | infinite-support `trim` behaviour |

**Tolerances tiered by estimator dependence:** `qi` uses type-7 quantiles → tight
`1e-8`; `hdi`/`mode`/`slab` depend on the KDE → looser, documented per fixture.
Pixel-level visual regression is rejected (brittle across Makie/font/backend
versions).

## Compatibility

- Julia `1.10+` (KernelDensity's floor; RandomDraws targets `1.9`, but this
  package's KDE dep raises the floor).
- Package extensions require Julia `1.9+` (satisfied).
- Pinned reference tooling: ggdist `3.3.3`, posterior `1.6.1`.

## Open questions / deferred

- Exact `rainclouds` attribute spellings — confirm against installed Makie
  during implementation.
- Discrete analytic distribution slabs — minimal/error-clear in v1; fuller
  support deferred.
- Custom AoG analyses — stretch goal beyond v1.
- Per-chain plotting — pooled in v1; per-chain deferred.
