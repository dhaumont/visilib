/*
Visilib, an open source library for exact visibility computation.
Copyright(C) 2021 by Denis Haumont

This file is part of Visilib.

Visilib is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Visilib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Visilib. If not, see <http://www.gnu.org/licenses/>
*/

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cstdlib>

#include "../demo/demo_helper.h"
#include "test_definition.h"

#include "geometry_convex_polygon.h"
#include "geometry_occluder_set.h"
#include "math_arithmetic.h"
#include "math_geometry.h"
#include "silhouette_mesh_face.h"
#include "math_vector_3.h"

using namespace visilib;
using namespace visilibDemo;

/* Consistency test for the visibility solvers.

The scene, the source parameters and the expected answer of each configuration
are loaded from the json files of the test case folder, the same files the
demo's write config feature produces, so the test runs exactly what can be
replayed and inspected in the interface.

Each configuration is solved twice, with the reference solver
(EXACT_APERTURE_FINDER) and with the solver under construction
(EXACT_SEQUENTIAL_SOLVER), for every combination of global scaling
{1e-5, 1, 1e5} with the silhouette optimisation and the early stop, on and
off. The scene scales together with its sources and the two flags only trade
speed for speed, so none of the three settings may change the answer: the
expectation of the case is asserted for each combination, and a single
deviating one fails the case. Each answer is compared to the expectation on
its own: the solvers are not compared to each other.

A sampling witness (probeWithSamples), third to the two solvers and sharing
no code with them, is kept around but sits idle: it uniformly samples
segments between the sources, intersected against the raw occluder triangles,
and can only ever prove visibility, since one unobstructed segment is a
certificate while the absence of one proves nothing. */

namespace
{
    const size_t cSampleCount = 50000;
    const double cBarycentricMargin = 1e-4;
    const double cBarycentricJitter = 1e-5;
    const double cSurfaceMargin = 1e-7;

    /* Column width of the result line: enough for FAILURE and for the variant counts. */
    const size_t cResultWidth = 8;

    /* One source primitive, ready for uniform sampling. Uses the library's
    GeometryConvexPolygon together with MathGeometry::uniformSampleTriangle for
    area-weighted sampling, while tracking barycentric coordinates for neighbourhood
    certification. Degenerate sources (point or segment) are handled as special cases. */

    struct SourceSampler
    {
        GeometryConvexPolygon mPolygon;
        std::vector<int> mCumulativeProbabilityTable;
        size_t mLookupSize;

        SourceSampler(const std::vector<float>& vertices)
            : mPolygon(&vertices[0], vertices.size() / 3), mLookupSize(2048)
        {
            if (mPolygon.getVertexCount() > 2)
            {
                std::vector<double> triangleFanAreas;
                MathGeometry::getTriangleFanAreas(mPolygon, triangleFanAreas);
                MathGeometry::computeCumulativeProbabilityLookupTable(mLookupSize, triangleFanAreas, mCumulativeProbabilityTable);
            }
        }

        size_t degreesOfFreedom() const
        {
            return mPolygon.getVertexCount() >= 3 ? 2 : mPolygon.getVertexCount() - 1;
        }
    };

    struct SourceSample
    {
        size_t mCorner[3];
        double mWeight[3];
        MathVector3d mPoint;
    };

    MathVector3d evaluate(const SourceSampler& sampler, const SourceSample& sample)
    {
        const std::vector<MathVector3d>& v = sampler.mPolygon.getVertices();

        MathVector3d result = v[sample.mCorner[0]];
        result *= sample.mWeight[0];

        for (size_t i = 1; i < 3; i++)
        {
            MathVector3d term = v[sample.mCorner[i]];
            term *= sample.mWeight[i];
            result += term;
        }

        return result;
    }

    /* Returns false when the draw is too close to the source border to be usable as
    the centre of a neighbourhood, in which case the caller must ignore it. */

    bool sampleSource(const SourceSampler& sampler, SourceSample& sample)
    {
        const size_t count = sampler.mPolygon.getVertexCount();
        if (count == 0)
        {
            return false;
        }

        const std::vector<MathVector3d>& v = sampler.mPolygon.getVertices();

        sample.mCorner[0] = 0;
        sample.mCorner[1] = 0;
        sample.mCorner[2] = 0;
        sample.mWeight[0] = 1.0;
        sample.mWeight[1] = 0.0;
        sample.mWeight[2] = 0.0;

        if (count == 1)
        {
            sample.mPoint = v[0];
            return true;
        }

        if (count == 2)
        {
            const double t = MathArithmetic<double>::getRandom();
            sample.mCorner[1] = 1;
            sample.mCorner[2] = 1;
            sample.mWeight[1] = t;
            sample.mWeight[0] = 1.0 - t;
            sample.mPoint = evaluate(sampler, sample);

            return std::fmin(sample.mWeight[0], sample.mWeight[1]) > cBarycentricMargin;
        }

        /* Select triangle from the fan using the library's cumulative probability table */
        size_t triangle = 0;
        if (!sampler.mCumulativeProbabilityTable.empty())
        {
            int randomIndex = (int)(MathArithmetic<float>::getRandom() * (sampler.mLookupSize - 1));
            triangle = (size_t)sampler.mCumulativeProbabilityTable[randomIndex];
        }

        /* Sample barycentric coordinates using library's uniformSampleTriangle */
        MathVector2d u(MathArithmetic<double>::getRandom(), MathArithmetic<double>::getRandom());
        MathVector3d bary = MathGeometry::uniformSampleTriangle(u);

        sample.mCorner[0] = 0;
        sample.mCorner[1] = triangle + 1;
        sample.mCorner[2] = triangle + 2;
        sample.mWeight[0] = bary.x;
        sample.mWeight[1] = bary.y;
        sample.mWeight[2] = bary.z;
        sample.mPoint = evaluate(sampler, sample);

        return bary.x > cBarycentricMargin && bary.y > cBarycentricMargin && bary.z > cBarycentricMargin;
    }

    /* The neighbourhood of a sample, all of it still stabbing the source: barycentric
    displacements for a polygon, along the segment for a segment. A point source has
    no such freedom, so plain 3D displacements are returned instead; they do not
    certify the existence of a stabbing line, they only check that the answer is
    transversal and not the result of a grazing miss. */

    std::vector<MathVector3d> neighbourhood(const SourceSampler& sampler, const SourceSample& sample, double distance)
    {
        std::vector<MathVector3d> points;
        points.push_back(sample.mPoint);

        const size_t freedom = sampler.degreesOfFreedom();

        if (freedom >= 1)
        {
            SourceSample offset = sample;

            for (size_t i = 1; i <= freedom; i++)
            {
                for (int sign = -1; sign <= 1; sign += 2)
                {
                    offset.mWeight[0] = sample.mWeight[0] - sign * cBarycentricJitter;
                    offset.mWeight[i] = sample.mWeight[i] + sign * cBarycentricJitter;

                    if (offset.mWeight[0] > 0.0 && offset.mWeight[i] > 0.0)
                    {
                        points.push_back(evaluate(sampler, offset));
                    }
                }

                offset.mWeight[0] = sample.mWeight[0];
                offset.mWeight[i] = sample.mWeight[i];
            }

            return points;
        }

        for (int axis = 0; axis < 3; axis++)
        {
            for (int sign = -1; sign <= 1; sign += 2)
            {
                points.push_back(MathVector3d(
                    sample.mPoint.x + (axis == 0 ? sign * distance : 0.0),
                    sample.mPoint.y + (axis == 1 ? sign * distance : 0.0),
                    sample.mPoint.z + (axis == 2 ? sign * distance : 0.0)));
            }
        }

        return points;
    }

    /* Moller-Trumbore, restricted to the open segment and to the interior of the
    triangle. Degenerate or parallel configurations report no hit, which is the
    conservative answer for a blocker: the neighbourhood requirement of the witness
    below is what stops such a sample from certifying visibility. */

    bool segmentHitsTriangle(const MathVector3d& origin, const MathVector3d& direction, double length,
                             const MathVector3d& v0, const MathVector3d& v1, const MathVector3d& v2)
    {
        const MathVector3d e1 = v1 - v0;
        const MathVector3d e2 = v2 - v0;

        const MathVector3d p = MathVector3d::cross(direction, e2);
        const double determinant = e1.dot(p);
        if (std::fabs(determinant) < 1e-30)
        {
            return false;
        }

        const MathVector3d s = origin - v0;
        const double u = s.dot(p) / determinant;
        if (u <= cBarycentricMargin || u >= 1.0 - cBarycentricMargin)
        {
            return false;
        }

        const MathVector3d q = MathVector3d::cross(s, e1);
        const double v = direction.dot(q) / determinant;
        if (v <= cBarycentricMargin || u + v >= 1.0 - cBarycentricMargin)
        {
            return false;
        }

        const double t = e2.dot(q) / determinant;
        const double gap = length * cSurfaceMargin;

        return t > gap && t < length - gap;
    }

    bool segmentHitsOccluders(GeometryOccluderSet* occluderSet, const MathVector3d& a, const MathVector3d& b)
    {
        MathVector3d direction = b - a;
        const double length = direction.getNorm();
        if (length <= 0.0)
        {
            return false;
        }
        direction *= (1.0 / length);

        for (size_t geometryId = 0; geometryId < occluderSet->getOccluderCount(); geometryId++)
        {
            std::vector<SilhouetteMeshFace>* faces = occluderSet->getOccluderConnectedFaces(geometryId);

            for (size_t faceIndex = 0; faceIndex < faces->size(); faceIndex++)
            {
                const SilhouetteMeshFace& face = (*faces)[faceIndex];

                for (size_t corner = 1; corner + 1 < face.getVertexCount(); corner++)
                {
                    if (segmentHitsTriangle(a, direction, length,
                                            convert<MathVector3d>(face.getVertex(0)),
                                            convert<MathVector3d>(face.getVertex(corner)),
                                            convert<MathVector3d>(face.getVertex(corner + 1))))
                    {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    struct WitnessResult
    {
        size_t mSampleCount;
        size_t mFreeSampleCount;
        bool mProvenVisible;
        MathVector3d mCertifierA;
        MathVector3d mCertifierB;
    };

    WitnessResult probeWithSamples(GeometryOccluderSet* occluderSet, const std::vector<float>& v0, const std::vector<float>& v1)
    {
        const SourceSampler source0(v0);
        const SourceSampler source1(v1);

        WitnessResult result = { cSampleCount, 0, false, MathVector3d(), MathVector3d() };

        for (size_t sample = 0; sample < cSampleCount; sample++)
        {
            SourceSample a;
            SourceSample b;
            if (!sampleSource(source0, a) || !sampleSource(source1, b))
            {
                continue;
            }

            const double length = (b.mPoint - a.mPoint).getNorm();
            if (length <= 0.0)
            {
                continue;
            }

            if (segmentHitsOccluders(occluderSet, a.mPoint, b.mPoint))
            {
                continue;
            }

            result.mFreeSampleCount++;

            const std::vector<MathVector3d> around0 = neighbourhood(source0, a, length * cSurfaceMargin);
            const std::vector<MathVector3d> around1 = neighbourhood(source1, b, length * cSurfaceMargin);

            bool neighbourhoodFree = true;
            for (size_t i = 0; i < around0.size() && neighbourhoodFree; i++)
            {
                for (size_t j = 0; j < around1.size() && neighbourhoodFree; j++)
                {
                    if ((around1[j] - around0[i]).getNorm() > 0.0 &&
                        segmentHitsOccluders(occluderSet, around0[i], around1[j]))
                    {
                        neighbourhoodFree = false;
                    }
                }
            }

            if (neighbourhoodFree)
            {
                result.mProvenVisible = true;
                result.mSampleCount = sample + 1;
                result.mCertifierA = a.mPoint;
                result.mCertifierB = b.mPoint;
                return result;
            }
        }

        return result;
    }

    /* ------------------------------------------------------------------- driver */

    VisibilityResult solveWith(GeometryOccluderSet* occluderSet, const std::vector<float>& v0, const std::vector<float>& v1,
                               const DemoConfiguration& configuration, VisibilityExactQueryConfiguration::SolverType solverType)
    {
        VisibilityExactQueryConfiguration query;
        query.detectApertureOnly = configuration.detectApertureOnly;
        query.hyperSphereNormalization = configuration.normalization;
        query.silhouetteOptimization = configuration.silhouetteOptimisation;
        query.precision = configuration.precisionType;
        query.tolerance = configuration.tolerance;
        query.confidenceValue = configuration.confidenceValue;
        query.solverType = solverType;

        /* The demo derives the library threshold from the source size and stores that
        radius in the case file; redo the same conversion here so the loaded file reproduces
        the minimumNormalizedApertureSize the case was pinned with. */

        const double side = 2.0 * configuration.scaling;
        const double aperture = configuration.minimumApertureSize;
        query.minimumNormalizedApertureSize = (aperture * aperture) / (side * side);

        return visilib::areVisible(occluderSet, &v0[0], v0.size() / 3, &v1[0], v1.size() / 3, query, nullptr);
    }

    const char* nameOf(VisibilityResult result)
    {
        switch (result)
        {
            case VISIBLE: return "VISIBLE";
            case HIDDEN:  return "HIDDEN";
            case UNKNOWN: return "UNKNOWN";
            default:      return "FAILURE";
        }
    }

    /* The bounds of the four loops runCase sweeps every case through: the arithmetic,
    the global scaling of the scene, the silhouette optimisation and the early stop.
    Each suite builds its own, so it can widen or narrow its sweep without touching
    the others. The switches are held as values, not counts: the suite sets both the
    order of the sweep and what on and off mean. The answer is a property of the
    geometry, and the scene scales with its sources, so none of these may change
    it. Neither may the arithmetic: the precision only chooses how the same scene is
    computed. */

    struct SweepBounds
    {
        std::vector<VisibilityExactQueryConfiguration::PrecisionType> mPrecisions;
        std::vector<float> mGlobalScalings;
        std::vector<bool> mSilhouetteOptimisations;
        std::vector<bool> mDetectApertureOnly;

        /* The number of scaling/optimisation/early-stop combinations, one arithmetic
        at a time, and of all combinations including the arithmetic. */

        size_t settingsCount() const
        {
            return mGlobalScalings.size() * mSilhouetteOptimisations.size() * mDetectApertureOnly.size();
        }

        size_t variantCount() const
        {
            return settingsCount() * mPrecisions.size();
        }
    };

    /* Every arithmetic this build offers, named as the case files name it. The answer
    cannot depend on how the scene is computed, so each of them must reproduce the pin:
    where a floating point answer diverges and an exact one does not, the failure is
    conditioning rather than logic. */

    std::vector<VisibilityExactQueryConfiguration::PrecisionType> allPrecisions()
    {
        return
        {
            VisibilityExactQueryConfiguration::DOUBLE,
            VisibilityExactQueryConfiguration::FLOAT
        #ifdef ENABLE_LEDA
            , VisibilityExactQueryConfiguration::LEDA_REAL
        #endif
        #ifdef ENABLE_GMP
            , VisibilityExactQueryConfiguration::GMP_FLOAT
            , VisibilityExactQueryConfiguration::GMP_RATIONAL
        #endif
        #ifdef ENABLE_MPFR
            , VisibilityExactQueryConfiguration::MPFR
        #endif
        };
    }

    const char* nameOf(VisibilityExactQueryConfiguration::PrecisionType precision)
    {
        switch (precision)
        {
            case VisibilityExactQueryConfiguration::FLOAT: return "flt";
            case VisibilityExactQueryConfiguration::DOUBLE: return "dbl";
        #ifdef ENABLE_LEDA
            case VisibilityExactQueryConfiguration::LEDA_REAL: return "leda";
        #endif
        #ifdef ENABLE_GMP
            case VisibilityExactQueryConfiguration::GMP_FLOAT: return "gmpf";
            case VisibilityExactQueryConfiguration::GMP_RATIONAL: return "gmpq";
        #endif
        #ifdef ENABLE_MPFR
            case VisibilityExactQueryConfiguration::MPFR: return "mpfr";
        #endif
            default: return "unknown";
        }
    }

    /* The two solvers every variant runs on: the reference (EXACT_APERTURE_FINDER)
    first, the solver under construction (EXACT_SEQUENTIAL_SOLVER) second. */
    const std::vector<VisibilityExactQueryConfiguration::SolverType> cSolvers =
    {
        VisibilityExactQueryConfiguration::EXACT_APERTURE_FINDER,
        VisibilityExactQueryConfiguration::EXACT_SEQUENTIAL_SOLVER
    };

    /* One label per solver, as the FAIL and the summary lines name them. */
    const char* cSolverNames[] = { "aperture", "sequential" };

    /* Prints 'label answer' for every solver, in the order of the cSolvers table. */

    void listAnswers(const std::vector<VisibilityResult>& answers, const char* lead)
    {
        std::cout << lead;

        for (size_t solverIndex = 0; solverIndex < cSolvers.size(); solverIndex++)
        {
            std::cout << (solverIndex == 0 ? " " : ", ") << cSolverNames[solverIndex] << " " << nameOf(answers[solverIndex]);
        }

        std::cout << std::endl;
    }

    /* Validates the answers of a variant against the expectation of the case: each
    solver's match is counted, per case and per precision, and a FAILURE or an answer
    beside the expectation is reported under the header. The answers are those of the
    cSolvers table, in order. Returns whether the variant failed, so the caller can
    count it and dump its configuration. */

    bool validate(const std::vector<VisibilityResult>& answers, VisibilityResult expected,
                  std::vector<size_t>& caseMatches, std::vector<size_t>& precisionMatches)
    {
        bool anyFailure = false;

        for (size_t solverIndex = 0; solverIndex < cSolvers.size(); solverIndex++)
        {
            if (answers[solverIndex] == expected)
            {
                caseMatches[solverIndex]++;
                precisionMatches[solverIndex]++;
            }

            if (answers[solverIndex] == FAILURE)
            {
                anyFailure = true;
            }
        }

        bool variantFailed = false;

        if (anyFailure)
        {
            listAnswers(answers, "    FAIL a solver answered FAILURE:");
            variantFailed = true;
        }

        if (expected == UNKNOWN)
        {
            listAnswers(answers, "    FAIL the case has no expectation,");
            variantFailed = true;
        }
        else
        {
            for (size_t solverIndex = 0; solverIndex < cSolvers.size(); solverIndex++)
            {
                if (answers[solverIndex] != expected)
                {
                    std::cout << "    FAIL " << cSolverNames[solverIndex] << " answered " << nameOf(answers[solverIndex])
                              << ", expected " << nameOf(expected) << std::endl;
                    variantFailed = true;
                }
            }
        }

        return variantFailed;
    }

    /* Runs one variant of the sweep of a case: builds the settings chosen by the
    indices, rebuilds the scene and the sources for them, solves the query with both
    solvers and compares each answer to the expectation of the case. The variant
    counter and the per-case and per-precision match counters are updated in place.
    Returns false, with the reason in errorMessage, when the scene could not be
    created and the case must stop. */

    bool test_solver(size_t& variant, size_t precisionIndex, size_t scalingIndex, size_t silhIndex, size_t earlyIndex,
                     const DemoConfiguration& configuration, const SweepBounds& sweep, VisibilityResult expected, const std::string& name, size_t caseIndex,
                     const char* tag, std::string& errorMessage, size_t& variantFailures, std::vector<size_t>& caseMatches, std::vector<size_t>& precisionMatches)
    {
        DemoConfiguration settings = configuration;
        settings.globalScaling = sweep.mGlobalScalings[scalingIndex];
        settings.silhouetteOptimisation = sweep.mSilhouetteOptimisations[silhIndex];
        settings.detectApertureOnly = sweep.mDetectApertureOnly[earlyIndex];
        settings.precisionType = sweep.mPrecisions[precisionIndex];

        /* Announced before the (possibly expensive) solves, so the line
        standing last in the output names the combination at work. */
        std::cout << "  [" << variant << "] gs=" << settings.globalScaling
                  << " silh=" << (settings.silhouetteOptimisation ? "on" : "off")
                  << " early=" << (settings.detectApertureOnly ? "on" : "off")
                  << " prec=" << nameOf(settings.precisionType) << std::endl;

        /* The global scaling enters both the scene and the sources: rebuild
        them for every combination. */
        std::vector<float> v0;
        std::vector<float> v1;

        DemoHelper::generatePolygon(v0, settings.vertexCount0, settings.scaling, settings.phi - 3.14519f, settings.globalScaling);
        DemoHelper::generatePolygon(v1, settings.vertexCount1, settings.scaling, settings.phi, settings.globalScaling);

        HelperTriangleMeshContainer* meshContainer = DemoHelper::createScene(settings.sceneIndex, settings.globalScaling);
        if (meshContainer == nullptr)
        {
            errorMessage = "could not create the scene of case '" + name + "'";
            return false;
        }

        GeometryOccluderSet* occluderSet = DemoHelper::createOccluderSet(meshContainer);

        /* The sequential solver draws from rand(), and the demo reseeds in
        initScene before every query. Do the same here so a variant answers
        the same thing whether it runs alone or after the others. */
        srand(0);

        /* Both solvers answer the same query: the reference first, the solver
        under construction second. The answers are then reported side by side. */
        std::vector<VisibilityResult> answers(cSolvers.size());
        for (size_t solverIndex = 0; solverIndex < cSolvers.size(); solverIndex++)
        {
            answers[solverIndex] = solveWith(occluderSet, v0, v1, settings, cSolvers[solverIndex]);
        }

        if (validate(answers, expected, caseMatches, precisionMatches))
        {
            variantFailures++;

            std::stringstream filename;
            filename << tag << "_" << caseIndex << "_" << variant << ".obj";
            DemoHelper::exportQueryToObj(filename.str(), v0, v1, *meshContainer, settings.globalScaling);
        }

        delete occluderSet;
        delete meshContainer;

        variant++;
        return true;
    }

    /* Returns false when the case could not be executed at all. The case announces
    itself before its sweep and closes with one result line, so a slow or crashing
    combination is always visible under its header. Failing variants are reported
    under the header as well, and the offending configuration is dumped as an .obj
    file named after the suite tag, for inspection. */

    bool runCase(const std::string& path, const std::string& name, size_t caseIndex, const SweepBounds& sweep, const char* tag, std::string& errorMessage, size_t& failureCount)
    {
        DemoConfiguration configuration;
        configuration.readConfig(path);

        /* readConfig keeps the defaults it was handed when the file is missing or
        unreadable. Those defaults are an empty source and an out of range scene, so a
        failed load stops right here rather than crashing in areVisible. */

        if (configuration.vertexCount0 == 0 || configuration.vertexCount1 == 0 || configuration.sceneIndex >= DemoConfiguration::COUNT)
        {
            errorMessage = "could not load the geometry of case '" + name + "' from " + path;
            return false;
        }

        /* Expected answer, witness checked by hand, pinned inside the case file under
        "expected": both solvers must reproduce it, in every combination of the sweep.
        It is what makes the test a regression net rather than a change detector: once
        the sequential solver is fixed it will be tempted to agree with itself, not
        with the truth. */

        const VisibilityResult expected = DemoConfiguration::resultFromName(configuration.expected);

        /* The header of the case, printed before any solving: whatever comes next,
        running or failing, belongs to it. */
        std::cout << std::endl;
        std::cout << "  " << name << std::endl;
        std::cout << "  scene=" << DemoConfiguration::sceneName(configuration.sceneIndex)
                  << " sources=" << configuration.vertexCount0 << "/" << configuration.vertexCount1
                  << " scaling=" << configuration.scaling
                  << " phi=" << configuration.phi
                  << " normalized=" << (configuration.normalization ? "on" : "off")
                  << " expected=" << nameOf(expected) << std::endl;

        std::vector<size_t> caseMatches(cSolvers.size(), 0);
        std::vector<size_t> precisionMatches(cSolvers.size(), 0);
        size_t variantFailures = 0;
        size_t variant = 0;

        for (size_t precisionIndex = 0; precisionIndex < sweep.mPrecisions.size(); precisionIndex++)
        {
            for (size_t scalingIndex = 0; scalingIndex < sweep.mGlobalScalings.size(); scalingIndex++)
            {
                for (size_t silhIndex = 0; silhIndex < sweep.mSilhouetteOptimisations.size(); silhIndex++)
                {
                    for (size_t earlyIndex = 0; earlyIndex < sweep.mDetectApertureOnly.size(); earlyIndex++)
                    {
                        if (!test_solver(variant, precisionIndex, scalingIndex, silhIndex, earlyIndex, configuration, sweep, expected, name, caseIndex,
                                         tag, errorMessage, variantFailures, caseMatches, precisionMatches))
                        {
                            return false;
                        }
                    }
                }
            }

            /* One arithmetic at a time: the sweep of the settings alone, so a precision
            that cannot hold the answer is recognisable without reading the variants. */
            std::cout << "    " << std::left << std::setw(5) << nameOf(sweep.mPrecisions[precisionIndex]);

            for (size_t solverIndex = 0; solverIndex < cSolvers.size(); solverIndex++)
            {
                std::cout << " " << cSolverNames[solverIndex] << " " << precisionMatches[solverIndex] << "/" << sweep.settingsCount();
                precisionMatches[solverIndex] = 0;
            }

            std::cout << std::endl;
        }

        /* The result line closing the case: the failing combinations have spoken for
        themselves under the header, the line gives what the sweep must show at a
        glance. The counts give, per solver, how many variants reproduced the
        expectation. */
        std::cout << "  RESULT = " << std::right << std::setw(cResultWidth) << nameOf(expected);

        for (size_t solverIndex = 0; solverIndex < cSolvers.size(); solverIndex++)
        {
            std::cout << " | " << std::setw(cResultWidth) << caseMatches[solverIndex] << "/" << sweep.variantCount();
        }

        std::cout << std::endl;

        failureCount += variantFailures;

        return true;
    }
}

bool SolverConsistencyTest(std::string& errorMessage)
{
    std::cout << "SolverConsistencyTest" << std::endl;

    /* The case files are the ones the demo itself writes, so any configuration seen in
    the interface can be dropped in the folder to join the sweep. */

    const std::string folder = DemoHelper::testCaseFolder();
    const std::vector<std::string> files = DemoHelper::listTestCaseFiles(folder);

    /* The consistency sweep: every arithmetic this build offers, the scene kept at
    its own scale, and each flag tried on then off. */

    const SweepBounds sweep = { allPrecisions(), { 1.0f }, { true, false }, { true, false } };

    size_t failureCount = 0;

    if (files.empty())
    {
        errorMessage = "no case file found in " + folder;
        std::cout << "    ERROR " << errorMessage << std::endl;
        return false;
    }

    for (size_t caseIndex = 0; caseIndex < files.size(); caseIndex++)
    {
        const std::string& name = files[caseIndex];

        if (!runCase(folder + "/" + name, name, caseIndex, sweep, "solverConsistency", errorMessage, failureCount))
        {
            std::cout << "    ERROR " << errorMessage << std::endl;
            return false;
        }
    }

    if (failureCount > 0)
    {
        std::cout << "SolverConsistencyTest FAILED: " << failureCount << " solver variants answered beside their expectation, over "
                  << files.size() << " cases and " << sweep.variantCount() << " combinations each" << std::endl;
        errorMessage = "solver inconsistency detected";
        return false;
    }

    std::cout << "SolverConsistencyTest SUCCESS" << std::endl;
    return true;
}

/* The same sweep run on the edge cases: the configurations parked in the
deactivate/on_the_edge subfolder, sources flush with an occluder, grazing
silhouettes and other near-degenerate settings, kept out of the main folder so
they do not slow the consistency sweep but still must reproduce their pin. */

bool SolverEdgeCases(std::string& errorMessage)
{
    std::cout << "SolverEdgeCases" << std::endl;

    /* The edge cases sit in a subfolder of the case folder, and the listing is
    not recursive, so the main sweep never sees them and this one sees nothing else. */

    const std::string folder = DemoHelper::testCaseFolder() + "/deactivate/on_the_edge";
    const std::vector<std::string> files = DemoHelper::listTestCaseFiles(folder);

    /* The edge sweep: same bounds as the consistency one, the edge is in the scenes,
    not in the settings. */

    const SweepBounds sweep = { allPrecisions(), { 1.0f }, { true, false }, { true, false } };

    size_t failureCount = 0;

    if (files.empty())
    {
        errorMessage = "no case file found in " + folder;
        std::cout << "    ERROR " << errorMessage << std::endl;
        return false;
    }

    for (size_t caseIndex = 0; caseIndex < files.size(); caseIndex++)
    {
        const std::string& name = files[caseIndex];

        if (!runCase(folder + "/" + name, name, caseIndex, sweep, "solverEdge", errorMessage, failureCount))
        {
            std::cout << "    ERROR " << errorMessage << std::endl;
            return false;
        }
    }

    if (failureCount > 0)
    {
        std::cout << "SolverEdgeCases FAILED: " << failureCount << " solver variants answered beside their expectation, over "
                  << files.size() << " cases and " << sweep.variantCount() << " combinations each" << std::endl;
        errorMessage = "solver inconsistency detected";
        return false;
    }

    std::cout << "SolverEdgeCases SUCCESS" << std::endl;
    return true;
}

/* The same sweep run on the precision cases: the configurations parked in the
precision_test subfolder, each one pinning a scene whose answer is meant to
question the arithmetic rather than the logic, so the per-precision lines of the
sweep are what the case is about. */

bool SolverPrecisionTests(std::string& errorMessage)
{
    std::cout << "SolverPrecisionTests" << std::endl;

    /* The precision cases sit in a subfolder of the case folder, and the listing is
    not recursive, so the main sweep never sees them and this one sees nothing else. */

    const std::string folder = DemoHelper::testCaseFolder() + "/precision_test";
    const std::vector<std::string> files = DemoHelper::listTestCaseFiles(folder);

    /* The precision sweep: same bounds as the consistency one for now, the cases
    will ask their own questions of the arithmetic. */

    const SweepBounds sweep = { allPrecisions(), { 1.0f }, { true, false }, { true, false } };

    size_t failureCount = 0;

    if (files.empty())
    {
        errorMessage = "no case file found in " + folder;
        std::cout << "    ERROR " << errorMessage << std::endl;
        return false;
    }

    for (size_t caseIndex = 0; caseIndex < files.size(); caseIndex++)
    {
        const std::string& name = files[caseIndex];

        if (!runCase(folder + "/" + name, name, caseIndex, sweep, "solverPrecision", errorMessage, failureCount))
        {
            std::cout << "    ERROR " << errorMessage << std::endl;
            return false;
        }
    }

    if (failureCount > 0)
    {
        std::cout << "SolverPrecisionTests FAILED: " << failureCount << " solver variants answered beside their expectation, over "
                  << files.size() << " cases and " << sweep.variantCount() << " combinations each" << std::endl;
        errorMessage = "solver inconsistency detected";
        return false;
    }

    std::cout << "SolverPrecisionTests SUCCESS" << std::endl;
    return true;
}
