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

#pragma once

#include <algorithm>
#include <limits>
#include <vector>

#include "geometry_convex_polygon.h"
#include "geometry_occluder_set.h"
#include "math_frame_3.h"
#include "math_geometry.h"
#include "math_vector_2.h"
#include "math_vector_3.h"
#include "silhouette_mesh_face.h"

namespace visilib
{
    /** @brief The receiver polygon standing for the whole receiver plane, when the shadows are traced on all of it.

    The queries need a receiver polygon: the lines they hold run from the light to it. To trace the shadows on the
    whole plane of the receiver, the polygon is enlarged to hold every shadow the occluders can cast on the plane:
    the convex hull, in the plane, of the given receiver and of the projections from the light of the corners of the
    bounding box of the occluders, widened by a margin. A corner the light does not project onto the plane (its ray
    goes away from it) casts a shadow without end: it is taken at a far point of its direction, cFarFactor times the
    size of the scene away, which is where the trace stops.*/
    class ShadowReceiver
    {
    public:
        static constexpr double cFarFactor = 10.0;   /**< @brief How far, in sizes of the scene, an unbounded shadow is followed*/
        static constexpr double cMargin = 0.01;      /**< @brief The relative widening of the enlarged receiver*/

        /** @brief Compute the enlarged receiver of the plane of someReceiverVertices, with the same orientation (the same plane
        normal). someLightVertices is the light; for a directional light, the projection is along aLightDirection instead.
        Return false when the receiver is not a polygon or the scene holds no occluder; the result is left empty then.*/
        static bool computeInfiniteReceiver(GeometryOccluderSet* aScene,
                                            const float* someLightVertices, size_t aLightVertexCount,
                                            const float* someReceiverVertices, size_t aReceiverVertexCount,
                                            bool aDirectional,
                                            std::vector<float>& someResult)
        {
            someResult.clear();
            if (aReceiverVertexCount < 3)
            {
                return false;
            }

            GeometryConvexPolygon receiver(someReceiverVertices, aReceiverVertexCount);
            GeometryConvexPolygon light(someLightVertices, aLightVertexCount);
            const MathPlane3d plane = receiver.getPlane();
            const MathVector3d receiverCenter = MathGeometry::getGravityCenter(receiver);
            const MathVector3d lightCenter = MathGeometry::getGravityCenter(light);
            const MathFrame3d frame(receiverCenter, plane.mNormal, receiver.getVertex(1) - receiver.getVertex(0));

            // the bounding box of the occluders
            MathVector3d boxMin(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
            MathVector3d boxMax = boxMin * -1.0;
            bool hasOccluder = false;
            for (size_t i = 0; i < aScene->getOccluderCount(); i++)
            {
                for (const SilhouetteMeshFace& face : *aScene->getOccluderConnectedFaces(i))
                {
                    for (size_t k = 0; k < face.getVertexCount(); k++)
                    {
                        const MathVector3d v = convert<MathVector3d>(face.getVertex(k));
                        boxMin = MathVector3d(std::min(boxMin.x, v.x), std::min(boxMin.y, v.y), std::min(boxMin.z, v.z));
                        boxMax = MathVector3d(std::max(boxMax.x, v.x), std::max(boxMax.y, v.y), std::max(boxMax.z, v.z));
                        hasOccluder = true;
                    }
                }
            }
            if (!hasOccluder)
            {
                return false;
            }

            // the size of the scene: the occluders, the light and the receiver
            double size = (boxMax - boxMin).getNorm();
            for (int i = 0; i < light.getVertexCount(); i++)
            {
                size = std::max(size, (light.getVertex(i) - receiverCenter).getNorm());
            }
            for (int i = 0; i < receiver.getVertexCount(); i++)
            {
                size = std::max(size, (receiver.getVertex(i) - receiverCenter).getNorm());
            }
            const double far = cFarFactor * size;

            std::vector<MathVector2d> points;
            for (int i = 0; i < receiver.getVertexCount(); i++)
            {
                points.push_back(frame.toLocal(receiver.getVertex(i)));
            }

            const MathVector3d direction = receiverCenter - lightCenter;
            for (int corner = 0; corner < 8; corner++)
            {
                const MathVector3d c((corner & 1) ? boxMax.x : boxMin.x, (corner & 2) ? boxMax.y : boxMin.y, (corner & 4) ? boxMax.z : boxMin.z);
                const int lightCount = aDirectional ? 1 : light.getVertexCount();
                for (int i = 0; i < lightCount; i++)
                {
                    // the ray of the light through the corner: from the light vertex, or along the direction of a directional light
                    const MathVector3d origin = aDirectional ? c - direction : light.getVertex(i);
                    const MathVector3d ray = c - origin;
                    points.push_back(projectRay(frame, origin, ray, far));
                }
            }

            std::vector<MathVector2d> hull = getConvexHull(points);
            if (hull.size() < 3)
            {
                return false;
            }

            // widen by the margin, around the center of the hull
            MathVector2d center = MathVector2d::Zero();
            for (const MathVector2d& p : hull)
            {
                center += p;
            }
            center *= 1.0 / double(hull.size());

            // the hull is counter clockwise in the frame, which is the orientation of the normal of the plane: keep the one of the receiver
            const bool reverse = !isCounterClockwise(receiver, frame);
            if (reverse)
            {
                std::reverse(hull.begin(), hull.end());
            }
            for (const MathVector2d& p : hull)
            {
                const MathVector3d v = frame.toWorld(center + (p - center) * (1.0 + cMargin));
                someResult.push_back(float(v.x));
                someResult.push_back(float(v.y));
                someResult.push_back(float(v.z));
            }
            return true;
        }

    private:
        /** @brief Where the ray from anOrigin along aRay meets the plane of aFrame, in its coordinates; the far point of its direction in the plane when it does not.*/
        static MathVector2d projectRay(const MathFrame3d& aFrame, const MathVector3d& anOrigin, const MathVector3d& aRay, double aFar)
        {
            const double denominator = aRay.dot(aFrame.getNormal());
            const double distance = aFrame.getSignedDistance(anOrigin);
            if (denominator != 0.0)
            {
                const double t = -distance / denominator;
                if (t > 0.0)
                {
                    const MathVector2d p = aFrame.toLocal(anOrigin + aRay * t);
                    if (p.getNorm() <= aFar)
                    {
                        return p;
                    }
                }
            }
            // the shadow goes on without end: follow the direction of the ray, projected in the plane
            MathVector2d direction(aRay.dot(aFrame.getAxisU()), aRay.dot(aFrame.getAxisV()));
            if (direction.getNorm() == 0.0)
            {
                return aFrame.toLocal(anOrigin);
            }
            direction.normalize();
            return aFrame.toLocal(anOrigin) + direction * aFar;
        }

        static bool isCounterClockwise(const GeometryConvexPolygon& aPolygon, const MathFrame3d& aFrame)
        {
            double area = 0.0;
            for (int i = 0; i < aPolygon.getVertexCount(); i++)
            {
                const MathVector2d a = aFrame.toLocal(aPolygon.getVertex(i));
                const MathVector2d b = aFrame.toLocal(aPolygon.getVertex((i + 1) % aPolygon.getVertexCount()));
                area += a.x * b.y - a.y * b.x;
            }
            return area > 0.0;
        }

        /** @brief The convex hull of the points, counter clockwise (Andrew's monotone chain).*/
        static std::vector<MathVector2d> getConvexHull(std::vector<MathVector2d> points)
        {
            std::sort(points.begin(), points.end(), [](const MathVector2d& a, const MathVector2d& b) { return a.x < b.x || (a.x == b.x && a.y < b.y); });
            auto cross = [](const MathVector2d& o, const MathVector2d& a, const MathVector2d& b) { return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x); };
            std::vector<MathVector2d> hull(2 * points.size());
            size_t k = 0;
            for (size_t i = 0; i < points.size(); i++)
            {
                while (k >= 2 && cross(hull[k - 2], hull[k - 1], points[i]) <= 0.0) k--;
                hull[k++] = points[i];
            }
            for (size_t i = points.size() - 1, t = k + 1; i > 0; i--)
            {
                while (k >= t && cross(hull[k - 2], hull[k - 1], points[i - 1]) <= 0.0) k--;
                hull[k++] = points[i - 1];
            }
            hull.resize(k > 0 ? k - 1 : 0);
            return hull;
        }
    };
}
