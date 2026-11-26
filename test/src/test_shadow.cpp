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
#include <iostream>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <map>
#include <tuple>
#include <algorithm>
#include <cmath>

#include "visilib.h"
#include "visibility_event.h"
#include "../demo/demo_helper.h"
#include "test_definition.h"

using namespace visilib;
using namespace visilibDemo;

namespace
{
    /* One light and receiver configuration, and what its discontinuities must be made of. */
    struct ShadowCase
    {
        DemoConfiguration::SCENE_TYPE mScene;
        size_t mLightVertices;
        float mLightSize;
        float mReceiverSize;
        bool mDirectional;
        bool mExpectEEE;        // at least one conic arc
        bool mOnlyLightVertex;  // a point light: every arc is a pencil of the light vertex
        bool mInfinite;         // traced on the whole receiver plane: some arc must leave the receiver polygon
    };

    const char* eventName(ShadowDiscontinuityArc::Event e)
    {
        switch (e)
        {
        case ShadowDiscontinuityArc::VE_LIGHT_VERTEX: return "VE_LIGHT_VERTEX";
        case ShadowDiscontinuityArc::VE_OCCLUDER_VERTEX: return "VE_OCCLUDER_VERTEX";
        case ShadowDiscontinuityArc::EEE: return "EEE";
        case ShadowDiscontinuityArc::DIRECTIONAL: return "DIRECTIONAL";
        }
        return "?";
    }

    /* The line L(t) of the family of an arc, as it was traced. */
    MathPlucker6d getArcLine(const ShadowDiscontinuityArc& arc, double t)
    {
        MathVector3d d = MathVector3d::Zero();
        MathVector3d m = MathVector3d::Zero();
        double power = 1.0;
        for (size_t k = 0; k < arc.mDirection.size(); k++)
        {
            d += arc.mDirection[k] * power;
            m += arc.mLocation[k] * power;
            power *= t;
        }
        return MathPlucker6d(d.x, d.y, d.z, m.x, m.y, m.z);
    }

    /* Check the arcs of one case. Each arc is sampled inside its range, and the line of each sample must:
    - meet every edge of its event (the Klein form vanishes),
    - cross the light and the receiver polygons (the critical line is a line from the light to the receiver),
    - have free lines arbitrarily close to it (it bounds a visible cell): one of the lines joining random points next
      to its two ends is not occluded, which the point to point query of the library tells. The free lines can be a
      thin wedge of the neighbourhood (a line through the corner of a hole, leaning on a farther edge, has a few
      percent of them), hence the number of tries.
    The trace must lie on the curve it is stored with: the line of a segment, the conic of an EEE event. */
    bool checkCase(const ShadowCase& c, std::string& message)
    {
        std::vector<float> v0, v1;
        DemoHelper::generatePolygon(v0, c.mLightVertices, c.mLightSize, -3.14519f, 1.0f);
        DemoHelper::generatePolygon(v1, 3, c.mReceiverSize, 0.0f, 1.0f);

        auto meshContainer = DemoHelper::createScene(c.mScene, 1.0f);
        GeometryOccluderSet* occluderSet = DemoHelper::createOccluderSet(meshContainer);

        VisibilityExactQueryConfiguration config;
        config.precision = VisibilityExactQueryConfiguration::DOUBLE;
        config.directionalLight = c.mDirectional;
        config.infiniteReceiver = c.mInfinite;

        HelperVisualDebugger debugger;
        debugger.getTimeout().setMaximumDuration(60.0);
        debugger.setDisplayStatistics(false);

        std::vector<ShadowDiscontinuityArc> arcs;
        const VisibilityResult result = computeShadowDiscontinuities(occluderSet, &v0[0], v0.size() / 3, &v1[0], v1.size() / 3, arcs, config, &debugger);

        std::ostringstream out;
        out << DemoConfiguration::sceneName(c.mScene) << " light=" << c.mLightVertices << (c.mDirectional ? " directional" : "")
            << " size=" << c.mLightSize << "/" << c.mReceiverSize << (c.mInfinite ? " infinite" : "") << ": ";

        bool ok = true;
        if (result == FAILURE || result == UNKNOWN)
        {
            out << "query failed";
            ok = false;
        }
        if (arcs.empty())
        {
            out << "no arc";
            ok = false;
        }

        // the reference reads the geometry the query used: snapped to its grid when the integer predicates were on
        for (float& c : v0) c = occluderSet->getGrid().snap(c);
        for (float& c : v1) c = occluderSet->getGrid().snap(c);
        GeometryConvexPolygon light(&v0[0], v0.size() / 3);
        GeometryConvexPolygon receiver(&v1[0], v1.size() / 3);
        const double tolerance = 1e-6 * c.mReceiverSize;
        VisibilityExactQueryConfiguration pointConfiguration;
        srand(1);

        size_t eee = 0;
        size_t outside = 0;
        for (const ShadowDiscontinuityArc& arc : arcs)
        {
            eee += arc.mEvent == ShadowDiscontinuityArc::EEE ? 1 : 0;
            if (c.mDirectional && arc.mEvent != ShadowDiscontinuityArc::DIRECTIONAL)
            {
                out << "a directional light gave a " << eventName(arc.mEvent) << " arc; ";
                ok = false;
            }
            if (c.mOnlyLightVertex && arc.mEvent != ShadowDiscontinuityArc::VE_LIGHT_VERTEX)
            {
                out << "a point light gave a " << eventName(arc.mEvent) << " arc; ";
                ok = false;
            }
            if ((arc.mType == ShadowDiscontinuityArc::CONIC) != (arc.mEvent == ShadowDiscontinuityArc::EEE))
            {
                out << "the trace of a " << eventName(arc.mEvent) << " arc has the wrong type; ";
                ok = false;
            }

            const std::vector<MathVector2d> polyline = arc.toPolyline(9);
            for (size_t i = 1; i + 1 < polyline.size(); i++)
            {
                const double t = arc.mT0 + (arc.mT1 - arc.mT0) * double(i) / 8.0;
                const MathVector2d& uv = polyline[i];
                const MathVector3d x = arc.toPoint3D(uv);
                const MathPlucker6d line = getArcLine(arc, t);
                const MathVector3d d = line.getDirection();

                // the stored curve
                if (arc.mType == ShadowDiscontinuityArc::SEGMENT)
                {
                    const MathVector3d& l = arc.mLineCoefficients;
                    const double distance = std::fabs(l.x * uv.x + l.y * uv.y + l.z) / std::sqrt(l.x * l.x + l.y * l.y);
                    if (distance > tolerance)
                    {
                        out << "a segment sample is off its line by " << distance << "; ";
                        ok = false;
                    }
                }
                else if (arc.mConic.getAlgebraicDistance(uv) > tolerance)
                {
                    out << "a conic sample is off its conic by " << arc.mConic.getAlgebraicDistance(uv) << "; ";
                    ok = false;
                }

                // the edges of the event
                for (size_t e = 0; e + 1 < arc.mEventEdges.size(); e += 2)
                {
                    const MathPlucker6d edge(arc.mEventEdges[e], arc.mEventEdges[e + 1]);
                    const double klein = std::fabs(line.dot(edge)) / (d.getNorm() * (arc.mEventEdges[e + 1] - arc.mEventEdges[e]).getNorm());
                    if (klein > 1e-7)
                    {
                        out << eventName(arc.mEvent) << " line misses an edge of its event by " << klein << "; ";
                        ok = false;
                    }
                }

                // the light and the receiver
                const MathPlane3d lightPlane = light.getPlane();
                const MathVector3d y = x - d * ((lightPlane.mNormal.dot(x) + lightPlane.d) / lightPlane.mNormal.dot(d));
                if (!c.mDirectional && c.mLightVertices >= 3 && !MathGeometry::isPointInsidePolygon(light, y, tolerance))
                {
                    out << eventName(arc.mEvent) << " line misses the light; ";
                    ok = false;
                }
                if (c.mInfinite)
                {
                    const MathPlane3d receiverPlane = receiver.getPlane();
                    if (std::fabs(receiverPlane.mNormal.dot(x) + receiverPlane.d) > tolerance)
                    {
                        out << eventName(arc.mEvent) << " trace is off the receiver plane; ";
                        ok = false;
                    }
                    outside += MathGeometry::isPointInsidePolygon(receiver, x, tolerance) ? 0 : 1;
                }
                else if (!MathGeometry::isPointInsidePolygon(receiver, x, tolerance))
                {
                    out << eventName(arc.mEvent) << " line misses the receiver; ";
                    ok = false;
                }

                // free lines next to it: the directional rays all have one direction, which a point query does not keep
                if (c.mDirectional)
                {
                    continue;
                }
                const MathVector3d from = c.mLightVertices == 1 ? light.getVertex(0) : y;
                bool free = false;
                for (int k = 0; k < 1024 && !free; k++)
                {
                    MathVector3d dy(rand() / (double)RAND_MAX - 0.5, rand() / (double)RAND_MAX - 0.5, rand() / (double)RAND_MAX - 0.5);
                    MathVector3d dx(rand() / (double)RAND_MAX - 0.5, rand() / (double)RAND_MAX - 0.5, rand() / (double)RAND_MAX - 0.5);
                    if (c.mLightVertices == 1)
                    {
                        dy = MathVector3d::Zero();
                    }
                    const MathVector3d a = from + dy * (2e-4 * c.mLightSize);
                    const MathVector3d b = x + dx * (2e-4 * c.mReceiverSize);
                    const float fa[3] = { (float)a.x, (float)a.y, (float)a.z };
                    const float fb[3] = { (float)b.x, (float)b.y, (float)b.z };
                    free = areVisible(occluderSet, fa, 1, fb, 1, pointConfiguration) == VISIBLE;
                }
                if (!free)
                {
                    out << eventName(arc.mEvent) << " line has no free line next to it; ";
                    ok = false;
                }
            }
        }
        if (c.mInfinite && outside == 0)
        {
            out << "no trace outside the receiver polygon; ";
            ok = false;
        }
        if (c.mExpectEEE && eee == 0)
        {
            out << "no EEE arc";
            ok = false;
        }

        out << " arcs=" << arcs.size() << " EEE=" << eee << (ok ? " OK" : " FAILED");
        std::cout << "  " << out.str() << std::endl;
        if (!ok)
        {
            message += out.str() + "\n";
        }

        return ok;
    }

    /* The complete set of discontinuities of a single planar occluder, against a reference built from its geometry.

    For an occluder lying in one plane, the critical lines of the shadow on the receiver are known in advance:
    - each edge of the outline of the occluder seen from each vertex of the light (VE, a vertex of the light),
    - each edge of the light seen through each corner of the outline (VE, a vertex of the occluder), a corner being an
      outline vertex where the outline turns.
    The edges inside the occluder, and the vertices where the outline goes on straight (the sides of a subdivided
    triangle are cut in collinear pieces), make none. The check goes both ways: the traces of the reference are covered
    by the arcs (nothing missing), and the arcs lie on the traces of the reference (nothing spurious). The samples close
    to the boundary of the receiver are left out: a trace along that boundary is a coincidence of the geometry. */
    bool checkPlanarOccluder(DemoConfiguration::SCENE_TYPE aScene, float aPhi, bool aSilhouetteOptimization, bool anIntegerPredicates,
                             VisibilityExactQueryConfiguration::PrecisionType aPrecision, std::string& message)
    {
        const float lightSize = 0.1f;
        const float receiverSize = 0.3f;
        std::vector<float> v0, v1;
        DemoHelper::generatePolygon(v0, 3, lightSize, aPhi - float(M_PI), 1.0f);
        DemoHelper::generatePolygon(v1, 3, receiverSize, aPhi, 1.0f);

        auto meshContainer = DemoHelper::createScene(aScene, 1.0f);
        GeometryOccluderSet* occluderSet = DemoHelper::createOccluderSet(meshContainer);

        VisibilityExactQueryConfiguration config;
        config.precision = aPrecision;
        config.silhouetteOptimization = aSilhouetteOptimization;
        config.integerPredicates = anIntegerPredicates;

        HelperVisualDebugger debugger;
        debugger.getTimeout().setMaximumDuration(60.0);
        debugger.setDisplayStatistics(false);

        std::vector<ShadowDiscontinuityArc> arcs;
        computeShadowDiscontinuities(occluderSet, &v0[0], 3, &v1[0], 3, arcs, config, &debugger);

        for (float& c : v0) c = occluderSet->getGrid().snap(c);
        for (float& c : v1) c = occluderSet->getGrid().snap(c);
        GeometryConvexPolygon light(&v0[0], 3);
        GeometryConvexPolygon receiver(&v1[0], 3);
        const MathPlane3d receiverPlane = receiver.getPlane();
        const MathVector3d receiverCenter = MathGeometry::getGravityCenter(receiver);

        // the outline: the edges of one face only; its corners: the outline vertices where two outline edges are not collinear
        typedef std::tuple<float, float, float> Key;
        std::map<std::pair<Key, Key>, int> edgeUse;
        for (size_t o = 0; o < occluderSet->getOccluderCount(); o++)
        {
            for (const SilhouetteMeshFace& face : *occluderSet->getOccluderConnectedFaces(o))
            {
                for (size_t k = 0; k < face.getVertexCount(); k++)
                {
                    const MathVector3f& a = face.getVertex(k);
                    const MathVector3f& b = face.getVertex((k + 1) % face.getVertexCount());
                    const Key ka(a.x, a.y, a.z), kb(b.x, b.y, b.z);
                    edgeUse[std::make_pair(std::min(ka, kb), std::max(ka, kb))]++;
                }
            }
        }
        std::vector<std::pair<MathVector3d, MathVector3d>> outline;
        for (const auto& e : edgeUse)
        {
            if (e.second == 1)
            {
                outline.push_back(std::make_pair(MathVector3d(std::get<0>(e.first.first), std::get<1>(e.first.first), std::get<2>(e.first.first)),
                                                 MathVector3d(std::get<0>(e.first.second), std::get<1>(e.first.second), std::get<2>(e.first.second))));
            }
        }
        std::vector<MathVector3d> corners;
        for (size_t i = 0; i < outline.size(); i++)
        {
            for (int side = 0; side < 2; side++)
            {
                const MathVector3d v = side ? outline[i].second : outline[i].first;
                const MathVector3d d1 = (side ? outline[i].first : outline[i].second) - v;
                for (size_t j = i + 1; j < outline.size(); j++)
                {
                    if (!(outline[j].first == v || outline[j].second == v))
                    {
                        continue;
                    }
                    const MathVector3d d2 = (outline[j].first == v ? outline[j].second : outline[j].first) - v;
                    if (MathVector3d::cross(d1, d2).getNorm() > 1e-6 * d1.getNorm() * d2.getNorm())
                    {
                        corners.push_back(v);
                    }
                }
            }
        }

        // the point of the receiver plane on the ray from a through b, beyond b
        auto project = [&](const MathVector3d& a, const MathVector3d& b, MathVector3d& x)
        {
            const MathVector3d d = b - a;
            const double denominator = receiverPlane.mNormal.dot(d);
            if (denominator == 0.0)
            {
                return false;
            }
            const double t = -(receiverPlane.mNormal.dot(a) + receiverPlane.d) / denominator;
            x = a + d * t;
            return t > 1.0;
        };
        // inside the receiver, away from its boundary
        auto insideReceiver = [&](const MathVector3d& x)
        {
            if (!MathGeometry::isPointInsidePolygon(receiver, x, 1e-9))
            {
                return false;
            }
            for (int k = 0; k < receiver.getVertexCount(); k++)
            {
                const MathVector3d a = receiver.getVertex(k);
                MathVector3d inward = MathVector3d::cross(receiverPlane.mNormal, receiver.getVertex((k + 1) % receiver.getVertexCount()) - a);
                inward.normalize();
                if (inward.dot(receiverCenter - a) < 0.0)
                {
                    inward *= -1.0;
                }
                if (inward.dot(x - a) < 1e-4)
                {
                    return false;
                }
            }
            return true;
        };
        // the distance, in the receiver plane, from x to the segment [a, b]
        auto distanceToSegment = [](const MathVector3d& x, const MathVector3d& a, const MathVector3d& b)
        {
            const MathVector3d e = b - a;
            const double l = e.getSquaredNorm();
            const double u = l > 0.0 ? std::max(0.0, std::min(1.0, (x - a).dot(e) / l)) : 0.0;
            return (x - (a + e * u)).getNorm();
        };

        // the arcs as polylines in 3D
        std::vector<std::vector<MathVector3d>> polylines;
        for (const ShadowDiscontinuityArc& arc : arcs)
        {
            std::vector<MathVector3d> polyline;
            for (const MathVector2d& p : arc.toPolyline(64))
            {
                polyline.push_back(arc.toPoint3D(p));
            }
            polylines.push_back(polyline);
        }
        auto covered = [&](const MathVector3d& x)
        {
            for (const auto& polyline : polylines)
            {
                for (size_t i = 1; i < polyline.size(); i++)
                {
                    if (distanceToSegment(x, polyline[i - 1], polyline[i]) < 1e-6)
                    {
                        return true;
                    }
                }
            }
            return false;
        };

        // the traces of the reference, as segments of the receiver plane
        std::vector<std::pair<MathVector3d, MathVector3d>> traces;
        size_t expected = 0;
        size_t missing = 0;
        for (int i = 0; i < light.getVertexCount(); i++)
        {
            for (const auto& e : outline)
            {
                MathVector3d a, b;
                if (project(light.getVertex(i), e.first, a) && project(light.getVertex(i), e.second, b))
                {
                    traces.push_back(std::make_pair(a, b));
                }
                for (int k = 1; k < 20; k++)
                {
                    MathVector3d x;
                    if (project(light.getVertex(i), e.first + (e.second - e.first) * (k / 20.0), x) && insideReceiver(x))
                    {
                        expected++;
                        missing += covered(x) ? 0 : 1;
                    }
                }
            }
        }
        for (const MathVector3d& corner : corners)
        {
            for (int i = 0; i < light.getVertexCount(); i++)
            {
                const MathVector3d l0 = light.getVertex(i);
                const MathVector3d l1 = light.getVertex((i + 1) % light.getVertexCount());
                MathVector3d a, b;
                if (project(l0, corner, a) && project(l1, corner, b))
                {
                    traces.push_back(std::make_pair(a, b));
                }
                for (int k = 1; k < 20; k++)
                {
                    MathVector3d x;
                    if (project(l0 + (l1 - l0) * (k / 20.0), corner, x) && insideReceiver(x))
                    {
                        expected++;
                        missing += covered(x) ? 0 : 1;
                    }
                }
            }
        }

        // the arcs on the reference
        size_t samples = 0;
        size_t spurious = 0;
        for (const auto& polyline : polylines)
        {
            for (size_t i = 1; i + 1 < polyline.size(); i += 8)
            {
                samples++;
                bool explained = false;
                for (const auto& trace : traces)
                {
                    if (distanceToSegment(polyline[i], trace.first, trace.second) < 1e-6)
                    {
                        explained = true;
                        break;
                    }
                }
                spurious += explained ? 0 : 1;
            }
        }

        const bool ok = expected > 0 && missing == 0 && spurious == 0;
        std::ostringstream out;
        out << DemoConfiguration::sceneName(aScene) << " phi=" << aPhi << " silhouette=" << (aSilhouetteOptimization ? "on" : "off")
            << (anIntegerPredicates ? " integer" : "") << (aPrecision == VisibilityExactQueryConfiguration::DOUBLE ? "" : " exact")
            << ": arcs=" << arcs.size() << " missing " << missing << "/" << expected << " spurious " << spurious << "/" << samples << (ok ? " OK" : " FAILED");
        std::cout << "  " << out.str() << std::endl;
        if (!ok)
        {
            message += out.str() + "\n";
        }
        return ok;
    }
}

bool ShadowDiscontinuityTest(std::string& message)
{
    const ShadowCase cases[] =
    {
        // scene                                     light  sizes        directional EEE    point  infinite
        { DemoConfiguration::SINGLE_TRIANGLE_01,     3,     0.12f, 0.12f, false,      false, false, false },
        { DemoConfiguration::SINGLE_TRIANGLE_02,     3,     0.12f, 0.12f, false,      false, false, false },
        { DemoConfiguration::SINGLE_TRIANGLE_03,     3,     0.12f, 0.12f, false,      false, false, false },
        { DemoConfiguration::SINGLE_TRIANGLE_04,     3,     0.12f, 0.12f, false,      false, false, false },
        { DemoConfiguration::SIMPLE_CUBE,            3,     1.0f,  1.0f,  false,      false, false, false },
        { DemoConfiguration::SLOT_OFF_AXIS_01,       3,     0.5f,  0.5f,  false,      true,  false, false },
        { DemoConfiguration::SLOT_OFF_AXIS_02,       3,     0.5f,  0.5f,  false,      true,  false, false },
        { DemoConfiguration::NOISY_CUBE,             3,     0.5f,  0.5f,  false,      true,  false, false },
        { DemoConfiguration::SINGLE_TRIANGLE_02,     1,     0.12f, 0.3f,  false,      false, true, false  },
        { DemoConfiguration::SLOT_OFF_AXIS_01,       1,     0.5f,  0.8f,  false,      false, true, false  },
        { DemoConfiguration::SINGLE_TRIANGLE_01,     3,     0.12f, 0.12f, true,       false, false, false },
        { DemoConfiguration::SLOT_OFF_AXIS_01,       3,     0.5f,  0.5f,  true,       false, false, false },
        { DemoConfiguration::SINGLE_TRIANGLE_01,     3,     0.12f, 0.05f, false,      false, false, true  },
        { DemoConfiguration::SLOT_OFF_AXIS_01,       3,     0.5f,  0.2f,  false,      true,  false, true  },
        { DemoConfiguration::SINGLE_TRIANGLE_02,     1,     0.12f, 0.05f, false,      false, true,  true  },
        { DemoConfiguration::SINGLE_TRIANGLE_01,     3,     0.12f, 0.05f, true,       false, false, true  },
    };

    bool ok = true;
    for (const ShadowCase& c : cases)
    {
        ok = checkCase(c, message) && ok;
    }

    const DemoConfiguration::SCENE_TYPE planarScenes[] =
    {
        DemoConfiguration::SINGLE_TRIANGLE_01, DemoConfiguration::SINGLE_TRIANGLE_02,
        DemoConfiguration::SINGLE_TRIANGLE_03, DemoConfiguration::SINGLE_TRIANGLE_04
    };
    for (DemoConfiguration::SCENE_TYPE scene : planarScenes)
    {
        for (float phi : { 0.3f, 1.1f, 1.9f })
        {
            for (bool silhouette : { false, true })
            {
                // the same reference with the geometry snapped to the integer grid and the exact predicates
                for (bool integer : { false, true })
                {
                    ok = checkPlanarOccluder(scene, phi, silhouette, integer, VisibilityExactQueryConfiguration::DOUBLE, message) && ok;
                }
#ifdef ENABLE_GMP
                // in GMP rational, every decision of the CSG is exact and its ties are broken by the simulation of simplicity
                ok = checkPlanarOccluder(scene, phi, silhouette, false, VisibilityExactQueryConfiguration::GMP_RATIONAL, message) && ok;
#endif
            }
        }
    }
    return ok;
}
