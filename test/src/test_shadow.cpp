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
    - have free lines arbitrarily close to it (it bounds a visible cell): one of a few lines joining points next to
      its two ends is not occluded, which the point to point query of the library tells.
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
                for (int k = 0; k < 64 && !free; k++)
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
    return ok;
}
