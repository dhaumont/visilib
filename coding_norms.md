# Visilib Coding Norms

Naming and style conventions deduced from the analysis of all project header files
(`visilib/*.h`, `visilib/*.hpp`, `demo/*.h`, `test/src/*.h`). External libraries
(`external/`, `demo/external/`) are excluded.

Each norm is stated as it is actually applied in the codebase, with source examples.

---

## 1. Files

| Element | Convention | Examples |
|---|---|---|
| Header name | `snake_case.h`, all lowercase | `math_vector_3.h`, `geometry_occluder_set.h` |
| Header name prefix | Domain category first: `math_`, `geometry_`, `plucker_`, `silhouette_`, `visibility_`, `helper_` | `plucker_polytope.h`, `helper_triangle_mesh.h` |
| Include guard | `#pragma once` (never `#ifndef` guards) | every header |
| Public umbrella header | `visilib.h` includes the API; `visilib.hpp` contains inline implementations included at its end | `visilib.h` → `#include "visilib.hpp"` |
| License header | GPL v3 comment block at top of every file | every header |

## 2. Namespace

- All library code lives in a single lowercase namespace: `namespace visilib { ... }`.
- The namespace body is fully indented (4 spaces).
- No nested sub-namespaces.

## 3. Classes and structs

| Element | Convention | Examples |
|---|---|---|
| Class name | `PascalCase`, prefixed with its domain category matching the file: `Math*`, `Geometry*`, `Plucker*`, `Silhouette*`, `Visibility*`, `Helper*` | `MathVector3_`, `GeometryAABB`, `PluckerPolytope`, `SilhouetteProcessor`, `VisibilityExactQuery`, `HelperStatisticCollector` |
| Concrete type-count suffix | Dimension appended to the geometric kind: `Vector2`, `Vector3`, `Plane3`, `Plucker6` | `MathVector3_`, `MathPlane3_`, `MathPlucker6` |
| Template class | Generic implementation classes end with a trailing underscore `_`; non-template classes do not | `MathVector3_`, `MathMatrix_`, `VisibilityExactQuery_` vs `PluckerPolytope`, `GeometryOccluderSet` |
| Template parameters | Single capital letters, semantically documented: `S` = scalar/coordinate representation, `P` = Plucker point type | `template<class S> class MathVector3_`, `template<class P, class S> class VisibilitySolver` |
| Interface vs. implementation | Abstract base class without suffix; templated implementation is the same name + `_` and inherits publicly | `VisibilityExactQuery` (interface, pure virtual) and `VisibilityExactQuery_ : public VisibilityExactQuery` |
| Nested helper type | Short `PascalCase` name, visibility where used | `struct Hash` inside `Silhouette` |
| `struct` vs `class` | `struct` for plain data / description records (public fields, may still have virtual methods); `class` for encapsulated types with private members | `struct GeometryDiscreteMeshDescription`, `struct VisibilityExactQueryConfiguration`, `struct VisibilityRay` vs `class GeometryOccluderSet` |
| Inheritance | Always `: public Base` | `class SilhouetteContainerEmbree : public SilhouetteContainer` |
| Virtual destructor | `virtual ~X() = default;` or empty inline body | `VisibilitySolver`, `GeometryOccluderSet` |
| Pure virtual method | Declared with `= 0` | `virtual VisibilityResult resolve() = 0;` |
| Override | `override` keyword on redefinitions | `getIndices(size_t aFace) const override` |
| Typedef of template | Concrete instantiations are `PascalCase` + lowercase scalar letter (`f`loat, `d`ouble, `i`nt) | `typedef MathVector3_<float> MathVector3f;`, `MathVector2i`, `MathMatrixd` |
| Typedef of external type | `PascalCase` alias with library prefix | `typedef leda::real MathLedaReal;` |

## 4. Class member variables

| Element | Convention | Examples |
|---|---|---|
| Private/protected members | `m` + `PascalCase` | `mVertices`, `mIndices`, `mQuery`, `mDebugger`, `mTimerFrequency`, `mConnectedFacesCache` |
| Bool members | `mIs*` / `mHas*` where the flag reads as a predicate | `mIsActive`, `mTimerIsRunning` |
| Collections | plural noun or `*Count`/`*s` | `mEdges`, `mSilhouettes`, `mCounts`, `mAvailableEdgeCount` |
| Public data members (structs) | plain `lowerCamelCase`, no `m` prefix | `vertexCount`, `faceCount`, `vertexArray`, `silhouetteOptimization`, `minimumNormalizedApertureSize` |
| Members mirroring a C API | C-style lowercase kept as-is | `org`, `dir`, `tnear`, `tfar` in `VisibilityRay` (Embree layout) |
| Public math components | short lowercase | `x, y, z` in `MathVector3_` |
| Member documentation | trailing Doxygen comment `/**< @brief ... */` aligned | `std::vector<MathVector3f> mVertices; /**< @brief The list of vertices */` |

## 5. Functions and methods

| Element | Convention | Examples |
|---|---|---|
| Name | `lowerCamelCase`, starts with a verb | `normalize()`, `addEdge()`, `extractStabbingLines()` |
| Getter | `get` + noun | `getNorm()`, `getMin()`, `getOccluderConnectedFaces()` |
| Count getter | `get<Collection>Count()` returning `size_t` | `getVertexCount()`, `getEdgeCount()`, `getExtremalStabbingLinesCount()` |
| Bool query | `is*` / `has*` / `contains*` (no `get`) | `isZero()`, `isFinite()`, `hasSomeEdgesCollapsed()`, `containsRealLines()` |
| Setter | `set*` | `setNeighbour()`, `setTimerIsRunning()`, `setGeometry()` |
| Computation | `compute*`, `find*`, `extract*`, `remove*`, `clear*`, `create*`, `build*`, `split*` | `computeExtremalStabbingLines()`, `findSceneIntersection()`, `removeDegeneratedTriangles()` |
| Free function | `lowerCamelCase` like methods | `areVisible()`, `arePolygonsVisible()`, `convert<To, From>()` |
| Static named value / factory | bare `PascalCase` (no verb), used for mathematical constants and tolerances | `MathVector3_::Zero()`, `X()`, `Y()`, `Z()`, `MathArithmetic<S>::Tolerance()` |
| Static utility class | all-static class named `Math*` / `Helper*` with only static methods | `MathArithmetic`, `MathPredicates`, `MathCombinatorial`, `MathGeometry` |
| Operator overloads | standard operators; stream insertion as `friend std::ostream& operator<<` inside the class | `operator+=`, `operator==`, `operator<<` in `MathVector3_` |
| Accessor pairs | const and non-const overloads returning container references | `const std::vector<...>& getVertices() const` and `std::vector<...>& getVertices()` |
| Inline definitions | small methods inline in the class; larger ones defined `inline` below the class in the same header or in `visilib.hpp` | `GeometryOccluderSet` methods, template bodies |

## 6. Function parameters

| Element | Convention | Examples |
|---|---|---|
| Name | `a` + `PascalCase` (article-style prefix) | `aFace`, `aMin`, `aPolyhedron`, `aConfiguration` |
| Vowel sound | `an` + `PascalCase` | `anEdgeIndex`, `anEnd`, `anEpsilon`, `anOffset` |
| Index/multi variants | suffix digits | `aPlane0`, `aPlane1`, `aFacetsDescription1`, `aV1`, `aV2` |
| Output parameters | regular names, passed by non-const reference | `MathVector3_<S>& aMin`, `std::vector<SilhouetteMeshFace>& aFaces` |
| Default values | allowed, typically `nullptr` or `0` | `const S& aDistance = 0`, `HelperVisualDebugger* debugger = nullptr` |
| `const` correctness | inputs are `const T&` or `const T*`; scalars by value | `isZero(S anEps)`, `addFace(const SilhouetteMeshFace& face)` |

## 7. Local variables

| Element | Convention | Examples |
|---|---|---|
| Primary style | `lowerCamelCase`, descriptive noun | `result`, `edgeHash`, `faceHash`, `frequency`, `geometryId`, `known`, `lines` |
| Secondary style | `my` + `PascalCase` (Hungarian-lite, in this codebase means "local temporary", used in inline/template code) | `myFaces`, `myBegin`, `myEnd`, `myEdgesMap`, `myFaceNumber`, `myTriangle` |
| Loop counters / iterators | single letters or plain names | `i`, `iter`, `line` |
| Short math temporaries | single letters matching math usage | `d` (norm), `inv`, `f`, `t`, `key` |
| Declaration | at first use, initialized | `MathVector3_ result(*this);` |

Both `lowerCamelCase` and `myX` are accepted; `lowerCamelCase` is preferred in newer/simple code, `myX` appears where many same-typed temporaries coexist.

## 8. Enumerations

| Element | Convention | Examples |
|---|---|---|
| Kind | unscoped `enum` (C-style), not `enum class` | `GeometryPositionType`, `VisibilityResult` |
| Type name | `PascalCase`, often suffixed `Type` for positional/categorical enums; result types named `...Result` | `TimerType`, `CounterType`, `PolytopePositionType`, `VisibilityResult` |
| Values | `UPPER_SNAKE_CASE` | `ON_BOUNDARY`, `ON_NEGATIVE_SIDE`, `VISIBILITY_QUERY`, `RAY_COUNT` |
| Explicit values | assigned when the numeric value matters | `ON_BOUNDARY = 1`, `BLOCKED = 0`, `UNDEFINED = -1` |
| Sentinel value | last entry `*_LAST` or `*_COUNT` to size arrays | `TIMER_LAST`, `COUNTER_LAST`, `PRECISION_COUNT`, `SOLVER_COUNT` |
| Nested config enums | declared inside the configuration struct that owns them | `VisibilityExactQueryConfiguration::PrecisionType`, `::SolverType` |
| Value documentation | trailing `/**< @brief ... */` | `ON_BOUNDARY = 1, /**< @brief The points lies on the geometric entity*/` |

Observed deviations (do not propagate): `HelperSyntheticMeshBuilder::SubdivisionPatternType` uses `PascalCase` values (`Regular`, `OneInside`); the demo uses an all-caps type name `DemoConfiguration::SCENE_TYPE`.

## 9. Macros and compile-time constants

| Element | Convention | Examples |
|---|---|---|
| Library macros | `UPPER_SNAKE_CASE` with `V_` prefix, defined inside `namespace visilib` | `V_ASSERT(X)`, `V_LOG`, `V_TRACE(EXPRESSION)` |
| Feature/optional macros | `ENABLE_<LIB>` uppercase | `ENABLE_LEDA`, `ENABLE_GMP`, `ENABLE_MPFR` |
| Diagnostic macros | descriptive uppercase | `INTEGRITY_CHECK`, `OUTPUT_DEBUG_FILE` |
| Numeric `#define` constants | `UPPER_SNAKE_CASE` | `HIT_LIST_LENGTH 16` |
| Null pointer | `nullptr` (never `NULL`/`0`) | throughout |

## 10. Documentation (Doxygen)

| Element | Convention | Examples |
|---|---|---|
| Class brief | `/** @brief ... */` block immediately above the declaration | every class |
| Method brief | `/** @brief ... */` above declaration; `@param name: description`, `@return: description` with the *bare* parameter name | `@param aVertex0: the first vertex of the edge` |
| Member/enum-value brief | trailing `/**< @brief ... */` | see sections 4 and 8 |
| Main page | `/**@mainpage ... */` in `visilib.h` | |

## 11. Formatting

- Indentation: 4 spaces, no tabs.
- Braces: Allman style (opening brace on its own line) for classes, functions, `if`, `for`, `while`; single-line bodies may stay on one line in math operator code (`MathVector3_& operator += (...) { x += v.x; ... }`).
- Member-init list after `:` with continuing lines indented:
  ```cpp
  VisibilitySolver(VisibilityExactQuery_<P, S>* aSolver)
      :mDebugger(nullptr),
      mQuery(aSolver)
  ```
- Access modifiers (`public:` / `protected:` / `private:`) flush with the class body indent.
- Pointer/reference stars attach to the type: `PluckerPolytope<P>* polytope`, `const MathVector3_& v`.
- Space before `()` of methods with no argument is often omitted: `getNorm()const` appears; prefer `getNorm() const`.
- `template<class S>` written without space between `class` and parameter in older code; both `template<class S>` and `template <class S>` occur — prefer the latter consistently.

## 12. Summary cheat-sheet

| Construct | Pattern | Example |
|---|---|---|
| File | `category_topic.h` | `math_vector_3.h` |
| Namespace | lowercase | `visilib` |
| Class (template impl.) | `CategoryThing[_]` | `MathVector3_` |
| Concrete typedef | `CategoryThingX` (x=f/d/i) | `MathVector3d` |
| Struct (data) | `CategoryThing` | `SilhouetteEdge` |
| Member | `mThing` | `mVertices` |
| Bool member | `mIsThing` / public `flagName` | `mIsActive`, `useEmbree` |
| Method | `verbNoun` / `getNoun` / `isNoun` | `computeExtremalStabbingLines()` |
| Static named value | `Name` | `Zero()`, `Tolerance()` |
| Parameter | `aThing` / `anThing` | `anEdgeIndex` |
| Local | `myThing` | `myEdgesMap` |
| Enum type | `ThingType` / `ThingResult` | `GeometryPositionType` |
| Enum value | `THING_VALUE` (+`_LAST`/`_COUNT` sentinel) | `TIMER_LAST` |
| Macro | `V_THING` / `ENABLE_THING` | `V_ASSERT` |
| Template param | `S`, `P`, `To`, `From` | `template<class P, class S>` |
