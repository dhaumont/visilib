# Visilib architecture

How a visibility query actually flows through the headers in `visilib/`.
Line references are to the current tree; they are the fastest way to check a
claim before acting on it.

## One `areVisible()` call

`visilib::areVisible()` is declared at `visilib/visilib.h:121` and *defined* in
`visilib/visilib.hpp` (`visilib.h` includes the `.hpp` at its bottom, so the
public header is never self-contained-compiled without the implementation).

```
areVisible(scene, v0, n0, v1, n1, config, debugger)          visilib.hpp:44
 ├─ validate (non-null verts; dynamic_cast scene — tautology)  visilib.hpp:47-62
 ├─ switch (config.precision) → new VisibilityExactQuery_<P,S> visilib.hpp:66-112
 │    FLOAT is the `default:` arm: an unknown enum silently becomes float
 ├─ query->attachVisualisationDebugger(debugger)               visilib.hpp:114
 ├─ query->arePolygonsVisible(...)                             visibility_exact_query.h:310
 └─ delete query                                               visilib.hpp:134
```

Pipeline inside `arePolygonsVisible` (`visilib/visibility_exact_query.h:310-385`),
each stage wrapped in a `HelperScopedTimer`:

| Stage | What happens |
|---|---|
| `createInitialPolygons` `:247-307` | builds two `GeometryConvexPolygon`s, synthesises a support plane for degenerate sources (point/segment), then mutually guard-band-clips them; clipping failure ⇒ `FAILURE` |
| point↔point fast path `:331-335` | plain ray test, no Plücker maths at all |
| `SilhouetteProcessor::init` + `extractAllSilhouettes` `:342-343` | convex hull ("shaft") of the two sources, then candidate silhouette edges per occluder |
| `mSilhouetteContainer->prepare()` `:347` | builds the ray-acceleration back-end |
| `PluckerPolytopeBuilder<P,S>::build` `:353` | the initial (minimal) polytope — the set of lines stabbing both sources; `appendPolytope` into the complex |
| `switch (config.solverType)` `:359-374` | **no `default:`** — a new enum value yields a null `solver` |
| `solver->resolve()` `:379` | see "Solvers" |

Ownership: everything is raw `new`/`delete` — there is not a single smart
pointer in `visilib/`, `demo/` or `test/`. `mScene` is a **borrowed**
`GeometryOccluderSet*` that the query never deletes
(`visibility_exact_query.h:192,237-244`), so the caller (and, transitively, the
mesh container behind it) must outlive the call. Shutdown order in the demo is
occluder set first, mesh container second, because the mesh descriptions hold
pointers into the container's vectors.

## Template parameters `P` and `S`

- `S` — the scalar back-end: `float`, `double`, `MathLedaReal`, `MathGmpFloat`,
  `MathGmpRational`, `MathMpfr`.
- `P` — the Plücker point type, always `MathPlucker6<S>`: two `MathVector3_<S>`
  (`mDirection` = d, `mLocation` = m = begin × end, `math_plucker_6.h:65-75`).
  The `<P,S>` pairing is hard-coded in `visilib.hpp`; nothing infers it.

A back-end must specialise `MathArithmetic<S>` (`math_arithmetic.h:37-77`) —
`Tolerance`, `getAbs`, `getSqrt`, `isFinite`, `to_double`, `bitsCount`,
`validDecimalsCount` are declaration-only in the primary template — and the type
itself needs the usual operator set plus `operator<<`
(`math_gmp_float.h` is the reference implementation). Adding one means five
edits: `math_<x>.h`, `MathArithmetic` specialisation, a `PrecisionType`
enumerator behind the matching `ENABLE_*` in `visilib.h`, a `case` in
`visilib.hpp`, and the `convert<>` / quadric-intersection specialisations in
`math_vector_3.h` / `math_geometry.h`.

Two scope limits worth knowing: only the Plücker-space CSG maths uses `S` —
scene vertices are always `float` and source polygons always
`MathVector3d` (`geometry_convex_polygon.h:92`) — and roughly thirty call sites
inside template code hard-code `MathArithmetic<double>` or `<float>` tolerances
(`math_geometry.h` throughout, `silhouette_processor.h:200,206`), so an "exact"
back-end is not exact everywhere.

## Data structures in Plücker space

- **`PluckerPolyhedron<P>`** (`plucker_polyhedron.h`) is *not* a closed solid: it
  is one grow-only arena that pools both Plücker points and splitting
  hyperplanes in a single `mLines` vector, plus cached per-point quadric side
  (`mQuadricRelativePositions`) and, for each point, the sorted list of
  hyperplane indices containing it (`mFacetsDescription`). `add()` returns the
  index and asserts that a normalized point really is normalized (`:141-150`).
- **`PluckerPolytope<P>`** (`plucker_polytope.h`) stores **only the skeleton**: a
  set of vertex *indices* into the polyhedron and a set of canonical `(min,max)`
  edges. Two vertices are adjacent iff their facet lists share ≥ 3 hyperplanes
  (`math_combinatorial.h:41-52`), and a split creates a vertex whose facets are
  `intersect(f1,f2) ∪ {new plane}` (`plucker_polytope_splitter.h:207-211`).
  `containsRealLines()` (`:489-507`) is the emptiness test that terminates the
  search: it is false when no edge crosses the Plücker quadric.
- **`PluckerPolytopeComplex<P>`** (`plucker_polytope_complex.h`) despite the name
  holds **no CSG tree**: one polyhedron + a flat vector of polytopes + a
  `mNeedRepack` flag. The dtor deletes every surviving pointer (`:113-120`), so
  `delete polytope` must always be paired with `setPolytope(i, nullptr)` and a
  later `repack()` — only `VisibilitySequentialSolver` plays that game. The CSG
  tree exists implicitly, either as the aperture finder's recursion (children
  are stack temporaries, deleted before returning) or as the sequential solver's
  evolving polytope list. `occlusionTreeNodeSymbol` strings (`"L"`/`"R"`
  suffixes) exist only for logging and are empty in a default build.

Sign decisions all funnel through
`MathPredicates::getRelativePosition(dot, epsilon)` → `GeometryPositionType`
(`ON_NEGATIVE_SIDE / ON_BOUNDARY / ON_POSITIVE_SIDE / ON_BOTH_SIDES`,
`geometry_position_type.h:26-33`). `ON_BOUNDARY` handling is where degeneracies
are resolved or silently ignored, so touch it carefully.

## Silhouettes

`SilhouetteProcessor::init` (`silhouette_processor.h:119-135`) builds the convex
hull of the two sources; extraction (`:219-301`) then flood-fills face adjacency
and keeps an edge only if it passes three cheap tests: inside the shaft
(conservative bounding-sphere test, `math_geometry.h:221-256`), the face lies
between the two source support planes, and — when `silhouetteOptimization` is on
— the edge is boundary or convex (`isPotentialSilhouetteEdge`, `:164-217`).
Result: one `Silhouette` per occluder with a list of `SilhouetteEdge`, each
lazily caching the index of its Plücker hyperplane in the shared polyhedron
(`silhouette.h:42-58`, filled in `visibility_aperture_finder.h:201-212`).

Non-obvious behaviour:

- With `silhouetteOptimization` **off**, the flood fill never propagates, so you
  get roughly one silhouette per face, many of them with zero edges.
- `GeometryConvexHullBuilder::build` returns `nullptr` on failure
  (`geometry_convex_hull.h:58-68`); the only alarm is a `V_ASSERT` at
  `silhouette_processor.h:130`, and a null hull makes every face fail the shaft
  test, so the query degenerates to `VISIBLE`.
- `SilhouetteProcessor` has no destructor, so `mConvexHull` leaks once per query.
- Ray tests go through `SilhouetteContainer` (`silhouette_container.h`), whose
  CPU path uses watertight `MathGeometry::hitsTriangle` / `hitsCylinder`.
  The `EMBREE` alternative (`silhouette_container_embree.h`) is guarded both at
  compile time and by `config.useEmbree` (`visibility_exact_query.h:224-233`),
  but `USE_EMBREE` is off, `FindEmbree.cmake` only probes a Windows path, and
  its `intersect(VisibilityRay*) override` no longer matches the base signature
  — treat that file as currently unbuildable.

## Solvers

`VisibilitySolver<P,S>` (`visibility_solver.h`) is a thin interface plus the
shared `extractStabbingLines` helper. The three implementations differ far more
than the enum names suggest:

| `SolverType` | Class | Reality |
|---|---|---|
| `EXACT_APERTURE_FINDER` (default) | `VisibilityApertureFinder` | exact DFS over the implicit CSG tree: `collectAllOccluders` casts one *representative line* through the current aperture (projected gravity centre of the polytope, `math_geometry.h:952-968`) ⇒ miss = `VISIBLE`; `!containsRealLines()` ⇒ `HIDDEN`; otherwise pick the best edge (`findNextEdge`), split, recurse. `detectApertureOnly` is what enables every early `return` and the full-skeleton path |
| `EXACT_SEQUENTIAL_SOLVER` | `VisibilitySequentialSolver` | batch: for each silhouette, classify every live polytope (`BLOCKED`/`OUTSIDE`/`INTERSECTS`) and split the intersecting ones. Ignores `detectApertureOnly`, never ray-casts, no early exit, ~25 unconditional `std::cout` traces per call, and prints nothing useful in CI |
| `MONTE_CARLO` | `VisibilityAggressiveSolver` | not a Plücker algorithm at all: area-and-cosine-weighted rejection sampling of segment pairs, `sampleCount = log(1-conf)/log(1-aperture)` (`visibility_aggressive_solver.h:113-119`). `mTolerance` is stored and never read, and `extractStabbingLines` is declared (`:43`) but never defined, so anything that starts calling it will fail to link |

Config fields by consumer: `hyperSphereNormalization` + `tolerance` +
`detectApertureOnly` → aperture finder; `hyperSphereNormalization` + `tolerance`
→ sequential; `detectApertureOnly` + `minimumNormalizedApertureSize` +
`confidenceValue` → aggressive. `silhouetteOptimization` is read by the query,
`useEmbree` by the container selection.

## Cross-cutting

- **Hypersphere normalization** (`MathPlucker6<S>::getNormalized`,
  `math_plucker_6.h:177-195`) reduces numerical error by keeping Plücker points
  on the unit hypersphere. It is applied at the *call sites* — polytope builder,
  splitter's new vertices, both solvers' hyperplanes, the representative line —
  and **not** inside `PluckerPolyhedron::add`, so any new insertion point must
  pre-normalize or the pooled coordinates become inconsistent. Note
  `createInitialPolygons`' `normalization` argument is unnamed/unused in the
  definition (`visibility_exact_query.h:188` vs `:247`).
- **Tolerance**: `getComputationTolerance<S>()` (`visilib.hpp:34-42`) uses
  `MathArithmetic<S>::Tolerance()` when `config.tolerance == -1` (the default
  sentinel) and `S(config.tolerance)` otherwise. Guard-band clipping always uses
  `MathArithmetic<double>` regardless of `S` (`visibility_exact_query.h:286,293`).
- **Debug switches**: `INTEGRITY_CHECK` *is* the definition of `V_ASSERT`
  (`visilib_core.h:45-49`) and additionally gates the O(V·F) `isValid()` calls —
  so disabling it removes assertions that some code paths treat as reachable
  behaviour (e.g. `plucker_polytope.h:502,613`). `FULL_VALIDITY_CHECK` changes
  one thing: duplicate-facet-vertex rejection in
  `PluckerPolyhedron::checkFacetsDescription` (`:218-221`). `OUTPUT_DEBUG_FILE`
  is pure `V_LOG` logging except that it is what makes child-tree symbols
  non-empty.
- **Statistics** (`helper_statistic_collector.h`): `TimerType` =
  `VISIBILITY_QUERY, RAY_INTERSECTION, SILHOUETTE_PROCESSING, POLYTOPE_SPLIT,
  POLYTOPE_BUILD, STABBING_LINE_EXTRACTION, OCCLUDER_TREATMENT`;
  `CounterType` = `RAY_COUNT, POLYTOPE_SPLIT_COUNT, OCCLUDER_TRIANGLE_COUNT`
  (the last counts shaft-intersecting faces, not occluders hit). Nested scopes
  on the same timer trip `V_ASSERT` (`:100-105`), so don't re-time a stage that
  is already timed higher up. `getTimingStats` returns `"Query"` in seconds and
  everything else as a percentage, which the demo then prints with the same
  `%.2f%%` format (`demo/demo_main.cpp:375-381`).
- **Visual debugger** (`helper_visual_debugger.h`) is a passive buffer, not a
  callback interface: five `add*/get*` pairs of `MathVector3f` (removed edges,
  stabbing lines, sampling lines, extremal stabbing lines, removed triangles)
  plus a recycled `GeometryTriangleMeshDescription`. Writers push pairs for
  "lines" and triples for triangles; `demo_debug_visualisation_gl.cpp` re-reads
  them with matching strides and `clear()` nulls the arrays (`:181-194`), so the
  buffers must never be read before the first write of a query.
- **Dead or misleading code** (verified, in case you go looking for it):
  AABB occluder culling is not implemented (`GeometryOccluderSet::mBoundingBoxes`,
  `hitsBoundingBox`, `mLastHit` have no readers/writers);
  `PluckerPolyhedron::mNormalizations` is written and never read;
  `Silhouette::Hash` is never instantiated; `PluckerPolytope` carries unused
  `mEdgesProcessed`, `mRadius`, `mSilhouettes`; `MathPlucker2` is unreachable
  from `areVisible()`; the depth-limit message says "200" for a limit of 2000.
