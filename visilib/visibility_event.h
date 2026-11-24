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

#include <vector>

#include "math_conic_2.h"
#include "math_frame_3.h"
#include "math_vector_2.h"
#include "math_vector_3.h"

namespace visilib
{
    /** @brief A 3D edge whose line is a hyperplane of the polyhedron: the feature a visibility event is made of.

    The edges of the light (query polygon 0), of the receiver (query polygon 1) and of the occluders are the lines the
    critical lines lean on. A degenerated query polygon (a point or a segment) is padded by the builder with synthetic
    edges, so that each vertex lies at two of them: their synthetic endpoints are not geometry, but the real vertex two
    of them share still is.*/
    struct VisibilityFeature
    {
        enum Kind
        {
            LIGHT_EDGE,     /**< @brief An edge of the light (query polygon 0)*/
            RECEIVER_EDGE,  /**< @brief An edge of the receiver (query polygon 1)*/
            OCCLUDER_EDGE,  /**< @brief A silhouette edge of an occluder*/
        };

        Kind mKind;
        MathVector3d mEndpoints[2];     /**< @brief The two endpoints of the edge*/
        bool mRealEndpoint[2];          /**< @brief Whether each endpoint is a vertex of the geometry, false for the padding of a degenerated polygon*/
        bool mSynthetic;                /**< @brief Whether the edge itself is padding: its points other than the real endpoints are not geometry*/

        /** @brief The index (0 or 1) of the real endpoint shared with another feature, -1 if none: the vertex of a VE event.*/
        int getSharedRealEndpoint(const VisibilityFeature& other) const
        {
            for (int i = 0; i < 2; i++)
            {
                if (!mRealEndpoint[i])
                {
                    continue;
                }
                for (int j = 0; j < 2; j++)
                {
                    if (other.mRealEndpoint[j] && mEndpoints[i] == other.mEndpoints[j])
                    {
                        return i;
                    }
                }
            }
            return -1;
        }

        bool isSource() const
        {
            return mKind == LIGHT_EDGE || mKind == RECEIVER_EDGE;
        }
    };

    /** @brief One discontinuity curve of the shadow on the receiver plane, traced by a free critical line family.

    The family is a set of lines with one degree of freedom, L(t) for t in [mT0, mT1], all of them unoccluded between the
    light and the receiver, and leaning on the edges of the event:
    - VE (a vertex and an edge): the lines through the vertex meeting the edge, a pencil in a plane. Its trace on the
      receiver plane is a segment.
    - EEE (three edges): the transversals of three skew lines, a regulus. Its trace on the receiver plane is a conic arc.
    - DIRECTIONAL (a directional light and an edge): the rays of the light meeting the edge, a strip of a plane. Its trace
      is a segment.

    The coordinates (u, v) of the trace are those of the receiver frame; the arc holds the polynomial coefficients of its
    family so that it can be sampled on the right branch of its conic (see toPolyline).*/
    struct ShadowDiscontinuityArc
    {
        enum Type
        {
            SEGMENT,  /**< @brief The trace is a segment: VE and DIRECTIONAL events*/
            CONIC     /**< @brief The trace is a conic arc: EEE events*/
        };

        enum Event
        {
            VE_LIGHT_VERTEX,    /**< @brief A vertex of the light and an occluder edge*/
            VE_OCCLUDER_VERTEX, /**< @brief A vertex of an occluder and an edge of the light or of an occluder*/
            EEE,                /**< @brief Three edges, at least two of them occluder edges*/
            DIRECTIONAL         /**< @brief The rays of a directional light meeting an occluder edge*/
        };

        Type mType;
        Event mEvent;

        MathFrame3d mFrame;                       /**< @brief The receiver frame the trace is written in*/
        MathVector3d mLineCoefficients;           /**< @brief SEGMENT: the trace lies on l . (u, v, 1) = 0*/
        MathConic2d mConic;                       /**< @brief CONIC: the trace lies on the conic*/
        MathVector2d mEndpoints[2];               /**< @brief The receiver points of the two extreme lines L(mT0) and L(mT1)*/

        std::vector<MathVector3d> mEventEdges;    /**< @brief The 3D edges of the event, two points each (the vertex of a VE event is an endpoint of two of them)*/

        double mT0;                               /**< @brief The first parameter of the family*/
        double mT1;                               /**< @brief The last parameter of the family*/
        std::vector<MathVector3d> mDirection;     /**< @brief The direction of L(t): sum of mDirection[k] t^k*/
        std::vector<MathVector3d> mLocation;      /**< @brief The location (moment) of L(t): sum of mLocation[k] t^k*/

        /** @brief The point where L(t) crosses the receiver plane, in receiver coordinates. False when L(t) is parallel to the plane.*/
        bool getPoint(double t, MathVector2d& aPoint) const
        {
            MathVector3d d = MathVector3d::Zero();
            MathVector3d m = MathVector3d::Zero();
            double power = 1.0;
            for (size_t k = 0; k < mDirection.size(); k++)
            {
                d += mDirection[k] * power;
                m += mLocation[k] * power;
                power *= t;
            }
            const double dn = d.dot(mFrame.getNormal());
            if (dn == 0.0)
            {
                return false;
            }
            // the point of the line closest to the origin, slid along the line to the plane
            const MathVector3d p = MathVector3d::cross(d, m) * (1.0 / d.getSquaredNorm());
            const MathVector3d x = p - d * (mFrame.getSignedDistance(p) / dn);
            aPoint = mFrame.toLocal(x);
            return true;
        }

        /** @brief The point of the receiver plane with the coordinates aPoint.*/
        MathVector3d toPoint3D(const MathVector2d& aPoint) const
        {
            return mFrame.toWorld(aPoint);
        }

        /** @brief Sample the trace with numSamples points (at least two), evenly in the parameter of the family.

        The samples are the receiver points of the lines of the family, which keeps them on the branch of the conic the
        family actually sweeps, and inside the parameter range: the implicit conic alone would not tell its branches
        apart. A line parallel to the receiver plane has no point and is skipped.*/
        std::vector<MathVector2d> toPolyline(size_t numSamples) const
        {
            std::vector<MathVector2d> result;
            const size_t n = numSamples < 2 ? 2 : numSamples;
            for (size_t i = 0; i < n; i++)
            {
                const double t = mT0 + (mT1 - mT0) * double(i) / double(n - 1);
                MathVector2d p;
                if (getPoint(t, p))
                {
                    result.push_back(p);
                }
            }
            return result;
        }
    };
}
