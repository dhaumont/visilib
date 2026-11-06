# AGENTS.md

Guidance for AI coding agents working in this repository.

## What this is

Visilib — exact from-region visibility computation in 3D (polygon/segment/point
queries) via CSG on polytopes in Plücker space. Header-only C++17 library
(`visilib/*.h*`, ~11k lines), plus an OpenGL/GLUT + Dear ImGui demo and a
hand-rolled test binary. Library code and demo code are separate CMake
subprojects; the only cross-dependency is that the test target compiles
`demo/demo_helper.cpp` as its scene fixture builder.

## Build (Linux)

`./generateProject.sh` is the real workflow — it configures + builds **both**
`debug/` (Debug) and `release/` (Release), and it is the only place that
provides `LEDA_DIR`:

```bash
./generateProject.sh
```

`configure.sh` (creates `build/`) and the README (`source configure.sh`,
`cd build`) are stale: they never set `LEDA_DIR`, so `find_package(LEDA REQUIRED)`
hard-fails in a fresh shell. Fix: `export LEDA_DIR=/path/to/LEDA-7` before
configuring, or change line 3 of `generateProject.sh`. `FindLEDA.cmake` looks
only at `$ENV{LEDA_DIR}` for `incl/LEDA/numbers/real.h` + `libleda.a`.

Other required system deps: CMake ≥ 3.24.2, freeglut, OpenGL, and GMP/MPFR
(`USE_GMP`/`USE_MPFR` default True and the demo/test targets link `gmp`/`mpfr`
unconditionally). CGAL is mentioned in the README but is not used anywhere in
the build or the sources.

Incremental builds (targets are `visilib`, `visilibDemo`, `visilibTest`):

```bash
cmake --build debug --target visilibTest -j8    # incremental, ~11 s after a header edit
export DISPLAY=:0 && debug/demo/visilibDemo
```

Build out of source only. Never `cmake .` at the repo root (that is what CI
does, and it litters the tree with `CMakeCache.txt`/`Makefiles`). `debug/`,
`release/`, `build/` are gitignored but hold ~115 MB of binaries — never
reference them as source.

## Compile-time switches

All in the top-level `CMakeLists.txt`, as plain `set()` (not cache options, so
`-DUSE_LEDA=OFF` on the command line has no effect — edit the file):
`USE_LEDA`, `USE_GMP`, `USE_MPFR` (default True) → `ENABLE_LEDA` /
`ENABLE_GMP` / `ENABLE_MPFR`; `USE_EMBREE` (default False) → `EMBREE` (the
Embree back-end in `silhouette_container_embree.h` is stale — its `intersect`
override no longer matches the base signature and `FindEmbree.cmake` only probes
a Windows path, so turning this on will not compile).
Also `INTEGRITY_CHECK` (on by default), and commented-out
`FULL_VALIDITY_CHECK` / `OUTPUT_DEBUG_FILE`. These are header-wide: the
`PrecisionType` enumerators in the public `visilib.h` are `#ifdef`-guarded, so
turning a flag off changes the public API, not just the build. Toggling any of
them recompiles the whole tree.

## Tests

`visilibTest` is the entire suite — one `main()` (`test/src/test_full.cpp`)
that runs `MathCombinatorialTest()` then `VisibilityTest()`. No test framework,
no registration macro, no per-case selection:

```bash
debug/test/visilibTest                                  # ~17 s, must exit 0
ctest --test-dir debug -R '^visilibTest$'
```

- Only the `visilibTest` CTest entry is real. CMake's `foreach` in
  `test/CMakeLists.txt` also registers `test_full`, `test_math`,
  `test_visibility`, `demo_helper`, which invoke a nonexistent `test_full`
  binary — a bare `ctest` always shows failures unrelated to your change.
- Run tests from a scratch directory: `VisibilityTest` writes ~50 `N_M.obj`
  debug dumps into the current working directory.
- To exercise a subset, edit the loop vectors (`precision`, `globalScaling`,
  `vertexCount`, `normalizations`, `silhouetteOptimization`) at the top of
  `VisibilityTest()` in `test/src/test_visibility.cpp` and rebuild. Revert
  before committing.
- New test: add `bool YourNewTest(std::string& errorMessage)` in `test/src/`,
  declare it in `test/src/test_definition.h`, call it from `main()` and
  `return 1` on false.
- `vertexCount` above 3 hits `displaySummary("Not supported with this library
  version"); break;`, which exits the switch but *not* the loop, so those cases
  silently re-run the config defaults. Don't "fix" it without checking intent.

## Layout / architecture

Flat directories with filename prefixes (no subdirectories): `math_*`
(vectors, Plücker types, predicates, `MathArithmetic<T>` traits) →
`geometry_*` (AABB, convex polygon/hull, meshes, occluder set) →
`silhouette_*` (silhouette generation, ray-casting back-ends) →
`plucker_*` (polytope/polyhedron pools, optimal builder, splitter) →
`visibility_*` (query orchestration and solver strategies) → `helper_*`
(stats, scene reading, synthetic mesh building, visual debugger) →
`visilib.h` / `visilib.hpp` (public surface).

Read `ARCHITECTURE.md` before changing anything in `visibility_*`, `plucker_*`
or `silhouette_*`: it documents the `areVisible()` stage order, the `MathArithmetic<S>`
back-end contract, the skeleton-only polytope representation, what each solver
actually reads from the configuration, and the ownership/leak traps. The short
version: `PluckerPolytopeComplex` is *not* a CSG tree (flat owning pointer list),
the CSG tree is the aperture finder's recursion, and everything is raw
`new`/`delete` with the scene pointer borrowed.

- `visilib/visilib.h` is the public API (`areVisible()`, `VisibilityResult`,
  `VisibilityExactQueryConfiguration`) and `#include`s `visilib.hpp` at its
  bottom, where the per-precision template instantiations live. Public headers
  carry Doxygen `/** @brief */` comments — keep them.
- The "library" is built as a `STATIC` target whose sources are all headers,
  with `LINKER_LANGUAGE CXX`. The `*_Src` lists in `visilib/CMakeLists.txt`
  exist for IDE `source_group` only: consumers compile through
  `include_directories(../visilib)` and `#include`, so an unlisted header still
  compiles (e.g. `visibility_sequential_solver.h` and
  `visibility_aggressive_solver.h` are absent from the lists). Add new headers
  to the right list anyway, and demo `.cpp` files to `Src` in
  `demo/CMakeLists.txt` (those really are compiled).
- Two dispatch points matter when extending: `configuration.precision` is
  switched in `visilib/visilib.hpp` (one `VisibilityExactQuery_<PluckerT,
  Scalar>` instantiation per arithmetic back-end), and
  `configuration.solverType` is switched at `visilib/visibility_exact_query.h:359`
  (no `default:` — a new enum value there is silently unsolved).
  A new numeric back-end means: `math_<type>.h`, a `MathArithmetic<T>`
  specialization in `math_arithmetic.h`, a `PrecisionType` enumerator guarded
  by the matching `ENABLE_*`, and a branch in `visilib.hpp`.
- Solver naming is misleading: `MONTE_CARLO` maps to
  `VisibilityAggressiveSolver`, `EXACT_APERTURE_FINDER` to
  `VisibilityApertureFinder`. `VisibilitySolver` is the abstract interface.
- The vendored ImGui and freeglut live under `demo/external/` — never edit
  them. Same for `visilib/external/convexHull/`.
- `doc/html/` and `doc/latex/` are committed Doxygen output. Regenerate with
  the `doxygen` binary against `doc/doxygen.cfg` (INPUT is `../visilib/` +
  `../README.md`, HTML/LATEX into `doc/`); otherwise don't touch.

## Conventions

- No formatter, linter, typecheck, pre-commit hook, or `.editorconfig` exists —
  don't invent one. Match the file you are editing: 4-space indent, allman-ish
  braces on types, `camelCase` members, `aParam`-style argument names,
  `mMember` for fields, tab/space mix inherited from the original author.
- Platform-specific code uses `#ifdef _WIN32` / `#ifdef __APPLE__`, with Linux
  as the `#else` branch (see the GLUT header includes in `demo/demo_main.cpp`).
- `test/`, `visilib/`, `demo/` and the root `CMakeLists.txt` open with the
  GPL-3.0 copyright header — reproduce it in new files.
- Untracked scratch files sit at the repo root (`a.out`, `build.txt`,
  `graphics.txt`, `imgui.ini`, `config_10.txt`, `patch.patch`,
  `revert.patch`, `visilib/silhouette_processor.h.back`). They are not
  gitignored and are not inputs — never stage them, and never add a `git add -A`.
- GPL-3.0 licensed; the demo loads `config.txt` from the current working
  directory and writes `config_N.txt` beside it.
