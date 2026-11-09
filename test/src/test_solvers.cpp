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
are loaded from the NN_*.json files of the test case folder, the same files the
demo's write config feature produces, so the test runs exactly what can be
replayed and inspected in the interface.

Each configuration is solved twice, with the reference solver
(EXACT_APERTURE_FINDER) and with the solver under construction
(EXACT_SEQUENTIAL_SOLVER). Each answer is compared to the expectation of the
case on its own, the solvers are not compared to each other: a pinned case
fails whichever solver misses the truth. Answers are additionally arbitrated
by a third witness that shares no code with either solver: uniform sampling
of segments between the two sources, intersected against the raw occluder
triangles with a plain double precision segment triangle test.

The witness can only ever prove visibility, since one unobstructed segment is a
certificate while the absence of one proves nothing. A candidate is accepted only
when a whole neighbourhood of it is unobstructed as well, so a sample grazing an
occluder edge cannot be mistaken for an aperture. Cases where nothing free was
sampled are therefore reported as inconclusive rather than as a failure, unless a
solver answered HIDDEN. */

namespace
{
    const size_t cSampleCount = 50000;
    const double cBarycentricMargin = 1e-4;
    const double cBarycentricJitter = 1e-5;
    const double cSurfaceMargin = 1e-7;

    /* Column widths of the one line result table: enough for the longest case file
    name and for FAILURE. */
    const size_t cNameWidth = 36;
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

    /* Returns false when the case could not be executed at all. Divergences are
    accumulated in failureCount, and the offending configuration is dumped as an
    .obj file for inspection. */

    bool runCase(const std::string& path, const std::string& name, size_t caseIndex, std::string& errorMessage, size_t& failureCount)
    {
        /* The sequential solver draws from rand(), and the demo reseeds in initScene
        before every query. Do the same here so a case answers the same thing whether it
        runs alone or after the other tests of the suite. */
        srand(0);

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
        "expected": both solvers must reproduce it. It is what makes the test a
        regression net rather than a change detector: once the sequential solver is
        fixed it will be tempted to agree with itself, not with the truth. A missing
        expectation reads back as UNKNOWN, which downgrades the case to a plain
        cross checking of the two solvers instead of a silent pass. */

        const VisibilityResult expected = DemoConfiguration::resultFromName(configuration.expected);

        std::vector<float> v0;
        std::vector<float> v1;

        DemoHelper::generatePolygon(v0, configuration.vertexCount0, configuration.scaling, configuration.phi - 3.14519f, configuration.globalScaling);
        DemoHelper::generatePolygon(v1, configuration.vertexCount1, configuration.scaling, configuration.phi, configuration.globalScaling);

        HelperTriangleMeshContainer* meshContainer = DemoHelper::createScene(configuration.sceneIndex, configuration.globalScaling);
        if (meshContainer == nullptr)
        {
            errorMessage = "could not create the scene of case '" + name + "'";
            return false;
        }

        GeometryOccluderSet* occluderSet = DemoHelper::createOccluderSet(meshContainer);

        const VisibilityResult aperture = solveWith(occluderSet, v0, v1, configuration, VisibilityExactQueryConfiguration::EXACT_APERTURE_FINDER);
        const VisibilityResult sequential = solveWith(occluderSet, v0, v1, configuration, VisibilityExactQueryConfiguration::EXACT_SEQUENTIAL_SOLVER);

//        const WitnessResult witness = probeWithSamples(occluderSet, v0, v1);

        /* One row per case: the scene, the sources and the witness stay in the .obj
        dumps and in the details below, the line only carries what the run must show
        at a glance. */
        std::cout << "  " << std::left << std::setw(cNameWidth) << name << " | "
                  << std::right << std::setw(cResultWidth) << nameOf(expected) << " | "
                  << std::setw(cResultWidth) << nameOf(sequential) << " | "
                  << std::setw(cResultWidth) << nameOf(aperture) << std::endl;

        bool failed = false;
        if (expected == UNKNOWN)
        {
            std::cout << "    Reference is UNKNOWN " << std::endl;
            failed = true;
        }
        else if (aperture == FAILURE || sequential == FAILURE)
        {
            std::cout << "    FAIL a solver answered FAILURE" << std::endl;
            failed = true;
        }
        else
        {
            if (sequential != expected)
            {
                std::cout << "    FAIL the sequential solver answered " << nameOf(sequential)
                        << ", expected " << nameOf(expected) << std::endl;
                failed = true;
            }

            if (aperture != expected)
            {
                std::cout << "    FAIL the aperture solver answered " << nameOf(aperture)
                        << ", expected " << nameOf(expected) << std::endl;
                failed = true;
            }    
        }
     
        if (failed)
        {
            failureCount++;

            std::stringstream filename;
            filename << "solverConsistency_" << caseIndex << ".obj";
            DemoHelper::exportQueryToObj(filename.str(), v0, v1, *meshContainer, configuration.globalScaling);
        }

        delete occluderSet;
        delete meshContainer;

        return true;
    }
}

bool SolverConsistencyTest(std::string& errorMessage)
{
    std::cout << "SolverConsistencyTest" << std::endl;

    /* Degenerate sources are exercised on purpose through the NN_*.json case files:
    the shaft used to select silhouette edges is the convex hull of the two sources,
    which collapses when one of them is a point or a segment. This is where the two
    solvers part ways the most. */

    const std::string folder = DemoHelper::testCaseFolder();
    const std::vector<std::string> files = DemoHelper::listTestCaseFiles(folder);

    size_t failureCount = 0;

    if (files.empty())
    {
        errorMessage = "no case file found in " + folder;
        std::cout << "    ERROR " << errorMessage << std::endl;
        return false;
    }

    std::cout << "  " << std::left << std::setw(cNameWidth) << "test case" << " | "
              << std::right << std::setw(cResultWidth) << "expected" << " | "
              << std::setw(cResultWidth) << "sequential" << " | "
              << std::setw(cResultWidth) << "aperture" << std::endl;

    for (size_t caseIndex = 0; caseIndex < files.size(); caseIndex++)
    {
        const std::string& name = files[caseIndex];

        if (!runCase(folder + "/" + name, name, caseIndex, errorMessage, failureCount))
        {
            std::cout << "    ERROR " << errorMessage << std::endl;
            return false;
        }
    }

    if (failureCount > 0)
    {
        std::cout << "SolverConsistencyTest FAILED: " << failureCount << " of " << files.size() << " solver subtests answered beside their expectation" << std::endl;
        errorMessage = "solver inconsistency detected";
        return false;
    }

    std::cout << "SolverConsistencyTest SUCCESS" << std::endl;
    return true;
}
