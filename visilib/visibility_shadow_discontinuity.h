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

#include <cstddef>
#include <vector>
#include <cmath>

#include "math_vector_2.h"
#include "math_vector_3.h"

namespace visilib
{
    /** @brief One arc of a shadow discontinuity, traced on the receiver polygon.

    The arc is a piece of the intersection between a shadow swath and the receiver plane. The swath
    is the surface swept by the one parameter family of stabbing lines carried by one edge of a
    visibility polytope:

    - a VE event (two of the three event lines meet at a vertex) sweeps a plane, and the trace is a
      SEGMENT;
    - an EEE event (three independent event lines) sweeps a regulus, and the trace is a CONIC.

    The curve is stored under implicit form in the coordinates (u, v) of the receiver plane basis
    (mOrigin, mE1, mE2), and bounded by the intersections of the two extremal stabbing lines
    supporting the polytope edge with the receiver plane (mP0, mP1).

    A CONIC also carries the rational parameterization of the family, used by pointAt and
    toPolyline: a member of the family is built from a point a(t) of the first event line as the
    line joining a(t) to c(t) = plane(a(t), l2) intersection l3. The arc is sampled in that
    parameter between the endpoints values mT0 and mT1, which keeps every sample on the arc itself:
    neither branch selection nor root solving on the implicit conic is asked of the caller. */

    struct ShadowDiscontinuityArc
    {
        /** @brief v * s, without disturbing the members: the vector product by a scalar of MathVector3_ is a non const member.*/
        static MathVector3d scale(MathVector3d v, double s)
        {
            v *= s;
            return v;
        }

        enum Type
        {
            SEGMENT,    /**< @brief Straight trace: the swath is a plane (VE event, or directional light)*/
            CONIC       /**< @brief Second order trace: the swath is a regulus (EEE event)*/
        };

        Type mType;

        MathVector2d mP0;       /**< @brief Endpoints of the arc, in receiver plane coordinates*/
        MathVector2d mP1;

        double mLine[3];        /**< @brief SEGMENT: implicit form l . (u, v, 1) = 0*/
        double mConic[6];       /**< @brief CONIC: implicit form A u^2 + B u v + C v^2 + D u + E v + F = 0*/

        size_t mFacets[3];      /**< @brief The polyhedron indices of the hyperplanes carrying the family*/
        size_t mFacetCount;
        bool mHasOccluderFacet; /**< @brief True when at least one facet is an occluder silhouette edge: only such arcs are shadow boundaries*/

        MathVector3d mOrigin;   /**< @brief Receiver plane basis: a point of the plane and two orthonormal directions. (mOrigin, mE1, mE2) is right handed: mE1 x mE2 is the plane normal*/
        MathVector3d mE1;
        MathVector3d mE2;

        MathVector3d mA[3];     /**< @brief The three event lines, as a point and a normalized direction each (the regulus rulings, order: parameterized line first)*/
        MathVector3d mU[3];
        double mT0;             /**< @brief Parameters of mP0 and mP1 on the family*/
        double mT1;

        /** @brief Project a 3D point of the receiver plane onto the (u, v) basis*/
        MathVector2d toPoint2D(const MathVector3d& aPoint3D) const
        {
            MathVector3d d = aPoint3D - mOrigin;
            return MathVector2d(d.dot(mE1), d.dot(mE2));
        }

        /** @brief Lift a 2D point of the receiver plane coordinates to 3D*/
        MathVector3d toPoint3D(const MathVector2d& aPoint2D) const
        {
            MathVector3d result(mOrigin);
            result += scale(mE1, aPoint2D.x);
            result += scale(mE2, aPoint2D.y);
            return result;
        }

        /** @brief Evaluate the implicit form of the curve at a point of the receiver plane coordinates*/
        double implicitValue(const MathVector2d& aPoint) const
        {
            if (mType == SEGMENT)
            {
                return mLine[0] * aPoint.x + mLine[1] * aPoint.y + mLine[2];
            }

            const double u = aPoint.x;
            const double v = aPoint.y;
            return mConic[0] * u * u + mConic[1] * u * v + mConic[2] * v * v + mConic[3] * u + mConic[4] * v + mConic[5];
        }

        /** @brief The point of the arc at the given parameter.

        SEGMENT: t in [0, 1] along the chord. CONIC: t is a value of the family parameter, and the point is the
        receiver crossing of the family member built at t. Returns a NaN point when the construction degenerates
        (the parameterized line and the second one become coplanar, or the member runs parallel to the receiver
        plane); the caller is asked to keep the previous sample in that case.*/
        MathVector2d pointAt(double t) const
        {
            if (mType == SEGMENT)
            {
                MathVector3d a(mP1.x - mP0.x, mP1.y - mP0.y, 0.0);
                MathVector3d p(mP0.x, mP0.y, 0.0);
                p += a * t;
                return MathVector2d(p.x, p.y);
            }

            MathVector3d a(mA[0]);
            a += scale(mU[0], t);

            // the plane spanned by the second line and a(t); its meeting point with the third line
            MathVector3d n = MathVector3d::cross(mU[1], a - mA[1]);
            double normal = n.normalize();
            if (normal < 1e-12)
            {
                return MathVector2d(std::nan(""), std::nan(""));
            }

            double denom = mU[2].dot(n);
            if (std::fabs(denom) < 1e-12)
            {
                return MathVector2d(std::nan(""), std::nan(""));
            }

            double s = (a - mA[2]).dot(n) / denom;
            MathVector3d c(mA[2]);
            c += scale(mU[2], s);

            // the family member a(t) - c(t), met with the receiver plane
            MathVector3d dir = c - a;
            MathVector3d normal3 = MathVector3d::cross(mE1, mE2);
            double dn = dir.dot(normal3);
            if (std::fabs(dn) < 1e-12)
            {
                return MathVector2d(std::nan(""), std::nan(""));
            }

            double tau = normal3.dot(mOrigin - a) / dn;
            MathVector3d q(a);
            q += dir * tau;
            return toPoint2D(q);
        }

        /** @brief The family parameter of the member of the CONIC family crossing the given receiver plane point.

        Used to clip the arc against the receiver polygon: a boundary crossing seen on the plane is turned back
        into a value of t, so the arc can be split in the parameter its sampling runs on. Returns a NaN when the
        construction degenerates; callers must test the value with std::isnan before using it.*/
        double eventParameterAt(const MathVector2d& aPoint) const
        {
            if (mType == SEGMENT)
            {
                MathVector3d chord(mP1.x - mP0.x, mP1.y - mP0.y, 0.0);
                double c2 = chord.dot(chord);
                if (c2 < 1e-24)
                {
                    return std::nan("");
                }
                MathVector3d d(aPoint.x - mP0.x, aPoint.y - mP0.y, 0.0);
                return d.dot(chord) / c2;
            }

            MathVector3d q = toPoint3D(aPoint);

            // the member through q: the plane (q, l2) met with l3 gives the second point c of it
            MathVector3d n = MathVector3d::cross(mU[1], q - mA[1]);
            double normal = n.normalize();
            double s = 0.0;
            MathVector3d c;
            if (normal >= 1e-12)
            {
                double denom = mU[2].dot(n);
                if (std::fabs(denom) >= 1e-12)
                {
                    s = (q - mA[2]).dot(n) / denom;
                    c = mA[2];
                    c += scale(mU[2], s);

                    MathVector3d dir = c - q;
                    double t = 0.0;
                    if (meetOnLine(mA[0], mU[0], q, dir, t))
                    {
                        return t;
                    }
                }
            }

            // the first choice grazed a degeneracy: swap the roles of the second and third lines
            n = MathVector3d::cross(mU[2], q - mA[2]);
            normal = n.normalize();
            if (normal < 1e-12)
            {
                return std::nan("");
            }

            double denom = mU[1].dot(n);
            if (std::fabs(denom) < 1e-12)
            {
                return std::nan("");
            }

            s = (q - mA[1]).dot(n) / denom;
            c = mA[1];
            c += scale(mU[1], s);

            MathVector3d dir = c - q;
            double t = 0.0;
            if (meetOnLine(mA[0], mU[0], q, dir, t))
            {
                return t;
            }

            return std::nan("");
        }

        /** @brief Sample the arc as a polyline of numSamples points in the receiver plane coordinates.

        The first and last points are exactly mP0 and mP1. Interior samples that hit a degenerate construction
        keep the previous point.*/
        std::vector<MathVector2d> toPolyline(size_t numSamples) const
        {
            std::vector<MathVector2d> result;

            if (numSamples < 2)
            {
                numSamples = 2;
            }

            result.push_back(mP0);
            MathVector2d last = mP0;

            for (size_t i = 1; i + 1 < numSamples; i++)
            {
                double t = mT0 + (mT1 - mT0) * (double)i / (double)(numSamples - 1);
                MathVector2d p = pointAt(t);

                if (std::isnan(p.x) || std::isnan(p.y))
                {
                    p = last;
                }

                result.push_back(p);
                last = p;
            }

            result.push_back(mP1);
            return result;
        }

        /** @brief The parameter of the meeting point of two lines assumed to intersect, read on the first one.

        Solves a(t) = x1 + t u1 = x2 + w dir in the least squares sense over the two unknowns; the residual is
        kept below tolerance since the lines are only supposed to meet, numerically. Returns false when the two
        are nearly parallel and ill conditioned.*/
        static bool meetOnLine(const MathVector3d& x1, const MathVector3d& u1, const MathVector3d& x2, const MathVector3d& dir, double& t)
        {
            MathVector3d r = x2 - x1;

            double e11 = u1.dot(u1);
            double e12 = -u1.dot(dir);
            double e22 = dir.dot(dir);

            double det = e11 * e22 - e12 * e12;
            if (e11 < 1e-24 || e22 < 1e-24 || std::fabs(det) < 1e-12 * e11 * e22)
            {
                return false;
            }

            double b1 = u1.dot(r);
            double b2 = -dir.dot(r);

            t = (e22 * b1 - e12 * b2) / det;
            return true;
        }
    };
}
