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
#include <cmath>
#include <vector>

#include "geometry_convex_polygon.h"
#include "math_arithmetic.h"
#include "math_combinatorial.h"
#include "math_geometry.h"
#include "math_plucker_6.h"
#include "math_plucker_traits.h"
#include "plucker_polyhedron.h"
#include "plucker_polytope.h"
#include "visibility_event_registry.h"
#include "visibility_shadow_discontinuity.h"

namespace visilib
{
    /** @brief Traces the shadow discontinuities a certified visible cell casts on the receiver polygon.

    The skeleton of the polytope carries the visibility events the trace is made of: an edge of the
    polytope joins two vertices sharing cCommonFacetsForEdge facets, and it is the one parameter
    family of stabbing lines lying on those hyperplanes. Where the family crosses the Plucker
    quadric it holds the real lines, and the two crossings are extremal stabbing lines whose
    receiver points bound one discontinuity arc; the family itself sweeps the shadow swath between
    them.

    The swath is decided by the event lines read off the shared facets:

    - three lines in general position (an EEE event): their transversals sweep a regulus. A point x
      lies on one of them when a direction d exists with ((x - a_i) x u_i) . d = 0 for the three
      lines, that is when the three rows are dependent:

          F(x) = det [ (x - a_i) x u_i ] = 0

      The cubic part of the determinant vanishes identically, so F is a quadratic form and its
      zero set is the quadric of the swath. Restricted to the receiver plane it gives the conic of
      the arc.
    - two of the lines meeting at a vertex (a VE event, source or occluder): the family is the
      pencil of that vertex over the remaining line, the swath is a plane, and the trace is a
      segment of that plane on the receiver.
    - the directional family (cCommonFacetsForEdge == 1): an edge of a screen cell is held by one
      carrier line and runs with the fixed direction of the light, so its swath is the plane of the
      two and its trace is a segment as well.

    Arcs whose facets hold no silhouette edge are bounds of the sources themselves on the receiver,
    not shadows, and are left out; so are the events built on the synthetic edges the builder pads a
    degenerated source with.

    Every value is computed in double precision, as the visualization path of the solver already
    does: the exact arithmetic decided which cells are visible, and what remains here is the shape
    of the boundary of those cells, read on the receiver plane. */

    template<class P, class S>
    class ShadowDiscontinuityTracer
    {
    public:
        /** @brief Trace the discontinuity arcs of one certified visible cell onto the receiver polygon.

        The arcs are appended to anOutput. The cell must be alive: the trace reads its skeleton, its
        facet descriptions and its crossing with the Plucker quadric in place.

        @param aPolyhedron the polyhedron holding the lines and the facet descriptions
        @param aPolytope the visible cell to trace
        @param aRegistry the provenance of the hyperplanes (may be incomplete: unknown facets classify geometrically)
        @param aSource the first query polygon (the emitter)
        @param aReceiver the second query polygon (the receiver): the trace is written in its plane
        @param aTolerance the tolerance of the computation
        @param anOutput the arcs accumulated so far, appended to */
        static void traceLeafPolytope(PluckerPolyhedron<P>* aPolyhedron,
                                      PluckerPolytope<P>* aPolytope,
                                      const VisibilityEventRegistry& aRegistry,
                                      const GeometryConvexPolygon& aSource,
                                      const GeometryConvexPolygon& aReceiver,
                                      S aTolerance,
                                      std::vector<ShadowDiscontinuityArc>& anOutput);

    private:
        /** @brief The receiver plane frame: a point, two orthonormal directions (e1 x e2 = normal) and the plane normal*/
        struct ReceiverFrame
        {
            MathVector3d mOrigin;
            MathVector3d mE1;
            MathVector3d mE2;
            MathVector3d mNormal;
            double mSize;     /**< @brief The greatest distance of a vertex of the polygon to the origin: the size the traces are read against*/
            double mEps;      /**< @brief A length tolerance scaled to the polygon size*/
        };

        // the homogeneous cross of two points of the plane coordinates: the line through them, l . (u, v, 1) = 0
        static MathVector3d homogeneousCross(const MathVector2d& a, const MathVector2d& b)
        {
            return MathVector3d(a.y - b.y, b.x - a.x, a.x * b.y - a.y * b.x);
        }

        // the value of the implicit arc form at a point of the plane coordinates
        static double implicitValueAt(const ShadowDiscontinuityArc& anArc, const MathVector2d& aPoint)
        {
            return anArc.implicitValue(aPoint);
        }

        // the roots of a x^2 + b x + c, kept inside [lo, hi]
        static void solveQuadratic(double a, double b, double c, double lo, double hi, std::vector<double>& roots)
        {
            double scale = std::max(std::fabs(a), std::max(std::fabs(b), std::fabs(c)));
            if (scale <= 0.0)
            {
                return;
            }

            const double eps = 1e-14;

            if (std::fabs(a) <= eps * scale)
            {
                if (std::fabs(b) <= eps * scale)
                {
                    return;
                }
                double r = -c / b;
                if (r >= lo && r <= hi)
                {
                    roots.push_back(r);
                }
                return;
            }

            double disc = b * b - 4.0 * a * c;
            if (disc < 0.0)
            {
                return;
            }

            disc = MathArithmetic<double>::getSqrt(disc);
            double r0 = (-b - disc) / (2.0 * a);
            double r1 = (-b + disc) / (2.0 * a);

            if (r0 >= lo && r0 <= hi)
            {
                roots.push_back(r0);
            }
            if (r1 >= lo && r1 <= hi)
            {
                roots.push_back(r1);
            }
        }

        // the meeting point of two lines which are known to cross, from a point and a direction each
        static bool intersectPoint(const MathVector3d& x1, const MathVector3d& u1, const MathVector3d& x2, const MathVector3d& u2, MathVector3d& out)
        {
            MathVector3d c = MathVector3d::cross(u1, u2);
            double c2 = c.dot(c);
            if (c2 <= 1e-24)
            {
                return false;
            }

            MathVector3d r = x2 - x1;
            MathVector3d rx = MathVector3d::cross(r, u2);
            double t = rx.dot(c) / c2;

            out = x1;
            out += ShadowDiscontinuityArc::scale(u1, t);
            return true;
        }

        // an axis least parallel to the normal, to build the frame of the receiver plane
        static MathVector3d leastParallelAxis(const MathVector3d& aNormal)
        {
            MathVector3d axes[3] = { MathVector3d::X(), MathVector3d::Y(), MathVector3d::Z() };
            int best = 0;
            double bestAbs = std::fabs(aNormal.dot(axes[0]));

            for (int i = 1; i < 3; i++)
            {
                double d = std::fabs(aNormal.dot(axes[i]));
                if (d < bestAbs)
                {
                    bestAbs = d;
                    best = i;
                }
            }

            return axes[best];
        }

        // frame of the receiver plane from its polygon: a point, an orthonormal pair, the normal,
        // and the mEps tolerance the 2D tests are scaled to
        static ReceiverFrame makeFrame(const GeometryConvexPolygon& aReceiver)
        {
            ReceiverFrame frame;
            frame.mSize = 0.0;
            frame.mEps = 0.0;

            MathVector3d n = aReceiver.getPlane().getNormal();
            double n2 = n.dot(n);
            if (n2 <= 0.0)
            {
                return frame;
            }

            frame.mNormal = ShadowDiscontinuityArc::scale(n, 1.0 / MathArithmetic<double>::getSqrt(n2));
            frame.mOrigin = ShadowDiscontinuityArc::scale(n, -aReceiver.getPlane().d / n2);

            MathVector3d axis = leastParallelAxis(frame.mNormal);
            frame.mE1 = MathVector3d::cross(frame.mNormal, axis);
            frame.mE1.normalize();
            frame.mE2 = MathVector3d::cross(frame.mNormal, frame.mE1);

            double size = 0.0;
            for (size_t i = 0; i < aReceiver.getVertexCount(); i++)
            {
                MathVector3d d = aReceiver.getVertex(i) - frame.mOrigin;
                size = std::max(size, d.getNorm());
            }

            frame.mSize = size;
            frame.mEps = size * 1e-9;
            return frame;
        }

        static MathVector2d project(const ReceiverFrame& frame, const MathVector3d& aPoint)
        {
            MathVector3d d = aPoint - frame.mOrigin;
            return MathVector2d(d.dot(frame.mE1), d.dot(frame.mE2));
        }

        static std::vector<MathVector2d> projectPolygon(const ReceiverFrame& frame, const GeometryConvexPolygon& aReceiver)
        {
            std::vector<MathVector2d> polygon;
            for (size_t i = 0; i < aReceiver.getVertexCount(); i++)
            {
                polygon.push_back(project(frame, aReceiver.getVertex(i)));
            }
            return polygon;
        }

        // membership of a point in the convex receiver polygon, boundary included up to the frame eps
        static bool insidePolygon(const std::vector<MathVector2d>& polygon, const ReceiverFrame& frame, const MathVector2d& aPoint)
        {
            const size_t count = polygon.size();
            if (count < 3)
            {
                return false;
            }

            double area2 = 0.0;
            for (size_t i = 0; i < count; i++)
            {
                const MathVector2d& a = polygon[i];
                const MathVector2d& b = polygon[(i + 1) % count];
                area2 += a.x * b.y - a.y * b.x;
            }

            const double sign = area2 >= 0.0 ? 1.0 : -1.0;
            const double eps = frame.mEps * (frame.mEps + 1.0);

            for (size_t i = 0; i < count; i++)
            {
                const MathVector2d& a = polygon[i];
                const MathVector2d& b = polygon[(i + 1) % count];

                double cross = (b.x - a.x) * (aPoint.y - a.y) - (b.y - a.y) * (aPoint.x - a.x);
                if (sign * cross < -eps)
                {
                    return false;
                }
            }

            return true;
        }

        /** @brief The conic of a regulus read on the receiver frame.

        A point x lies on a transversal of the three lines (u_i, m_i) when the rows (x x u_i) - m_i
        are dependent: F(x) = det[...] . The apparent cubic part of the determinant vanishes
        identically ((x x u2) x (x x u3) = x [x u2 u3] is carried by x, orthogonal to x x u1), what
        remains is the quadratic form of the quadric of the swath:

          F(x) = x.S x + w.x + c

        with the pairs (linear rows i, j, constant row k) taken with the parity of the permutation
        of (0, 1, 2): (0, 1|2) +, (0, 2|1) -, (1, 2|0) +, and the lone linear rows i read against
        the cross of the two constant rows of the other ones. Restricted to the frame
        (x = origin + u e1 + v e2), it gives the coefficients (A, B, C, D, E, F) of

          A u^2 + B u v + C v^2 + D u + E v + F = 0

        normalized by the greatest coefficient. */
        static bool buildConic(const MathVector3d u[3], const MathVector3d m[3], const ReceiverFrame& frame, double coeffs[6])
        {
            // A[i] is the matrix of x -> x x u_i ; b[i] = -m_i
            double Ai[3][3][3];
            double bi[3][3];

            for (int i = 0; i < 3; i++)
            {
                const MathVector3d& ui = u[i];
                Ai[i][0][0] = 0.0;        Ai[i][0][1] = ui.z;    Ai[i][0][2] = -ui.y;
                Ai[i][1][0] = -ui.z;      Ai[i][1][1] = 0.0;     Ai[i][1][2] = ui.x;
                Ai[i][2][0] = ui.y;       Ai[i][2][1] = -ui.x;   Ai[i][2][2] = 0.0;

                bi[i][0] = -m[i].x;
                bi[i][1] = -m[i].y;
                bi[i][2] = -m[i].z;
            }

            double Q[4][4];
            for (int r = 0; r < 4; r++)
            {
                for (int c = 0; c < 4; c++)
                {
                    Q[r][c] = 0.0;
                }
            }

            const int pairs[3][4] = { { 0, 1, 2, 1 }, { 0, 2, 1, -1 }, { 1, 2, 0, 1 } };

            for (int t = 0; t < 3; t++)
            {
                const int i = pairs[t][0];
                const int j = pairs[t][1];
                const int k = pairs[t][2];
                const double sign = (double)pairs[t][3];

                // skew of b_k : [b]x v = b x v
                double bkx[3][3] =
                {
                    { 0.0, -bi[k][2], bi[k][1] },
                    { bi[k][2], 0.0, -bi[k][0] },
                    { -bi[k][1], bi[k][0], 0.0 }
                };

                // T = A_i^T [b_k]x A_j  (row r of A_i^T is column r of A_i)
                double T[3][3];
                for (int r = 0; r < 3; r++)
                {
                    for (int c = 0; c < 3; c++)
                    {
                        double sum = 0.0;
                        for (int a = 0; a < 3; a++)
                        {
                            for (int bq = 0; bq < 3; bq++)
                            {
                                sum += Ai[i][a][r] * bkx[a][bq] * Ai[j][bq][c];
                            }
                        }
                        T[r][c] = -sign * sum;
                    }
                }

                for (int r = 0; r < 3; r++)
                {
                    for (int c = 0; c < 3; c++)
                    {
                        Q[r][c] += 0.5 * (T[r][c] + T[c][r]);
                    }
                }
            }

            // linear part: w = sum_i A_i^T (b_j x b_k), cyclic
            const int cyc[3][3] = { { 0, 1, 2 }, { 1, 2, 0 }, { 2, 0, 1 } };
            double w[3] = { 0.0, 0.0, 0.0 };

            for (int t = 0; t < 3; t++)
            {
                const int i = cyc[t][0];
                const int j = cyc[t][1];
                const int k = cyc[t][2];

                double cx = bi[j][1] * bi[k][2] - bi[j][2] * bi[k][1];
                double cy = bi[j][2] * bi[k][0] - bi[j][0] * bi[k][2];
                double cz = bi[j][0] * bi[k][1] - bi[j][1] * bi[k][0];

                for (int r = 0; r < 3; r++)
                {
                    w[r] += Ai[i][0][r] * cx + Ai[i][1][r] * cy + Ai[i][2][r] * cz;
                }
            }

            for (int r = 0; r < 3; r++)
            {
                Q[r][3] = 0.5 * w[r];
                Q[3][r] = 0.5 * w[r];
            }

            // constant part: b_0 . (b_1 x b_2)
            Q[3][3] = bi[0][0] * (bi[1][1] * bi[2][2] - bi[1][2] * bi[2][1])
                    - bi[0][1] * (bi[1][0] * bi[2][2] - bi[1][2] * bi[2][0])
                    + bi[0][2] * (bi[1][0] * bi[2][1] - bi[1][1] * bi[2][0]);

            // restrict: C3 = B^T Q B, the columns of B being (e1, 0), (e2, 0), (origin, 1)
            double B[4][3] =
            {
                { frame.mE1.x, frame.mE2.x, frame.mOrigin.x },
                { frame.mE1.y, frame.mE2.y, frame.mOrigin.y },
                { frame.mE1.z, frame.mE2.z, frame.mOrigin.z },
                { 0.0,         0.0,         1.0              }
            };

            double C3[3][3];
            for (int p = 0; p < 3; p++)
            {
                for (int q = 0; q < 3; q++)
                {
                    double sum = 0.0;
                    for (int r = 0; r < 4; r++)
                    {
                        for (int c = 0; c < 4; c++)
                        {
                            sum += B[r][p] * Q[r][c] * B[c][q];
                        }
                    }
                    C3[p][q] = sum;
                }
            }

            coeffs[0] = C3[0][0];
            coeffs[1] = 2.0 * C3[0][1];
            coeffs[2] = C3[1][1];
            coeffs[3] = 2.0 * C3[0][2];
            coeffs[4] = 2.0 * C3[1][2];
            coeffs[5] = C3[2][2];

            double maxCoeff = 0.0;
            for (int i = 0; i < 6; i++)
            {
                maxCoeff = std::max(maxCoeff, std::fabs(coeffs[i]));
            }

            // a family with no second order part left: the three lines do not bound a regulus
            if (!(maxCoeff > 1e-300))
            {
                return false;
            }

            for (int i = 0; i < 6; i++)
            {
                coeffs[i] /= maxCoeff;
            }

            return true;
        }

        // the value of the arc's implicit form and its directional derivative on the plane coordinates
        static double implicitAt(const ShadowDiscontinuityArc& arc, double u, double v)
        {
            if (arc.mType == ShadowDiscontinuityArc::SEGMENT)
            {
                return arc.mLine[0] * u + arc.mLine[1] * v + arc.mLine[2];
            }

            return arc.mConic[0] * u * u + arc.mConic[1] * u * v + arc.mConic[2] * v * v
                 + arc.mConic[3] * u + arc.mConic[4] * v + arc.mConic[5];
        }

        static double implicitDerivative(const ShadowDiscontinuityArc& arc, double u, double v, double eu, double ev)
        {
            if (arc.mType == ShadowDiscontinuityArc::SEGMENT)
            {
                return arc.mLine[0] * eu + arc.mLine[1] * ev;
            }

            const double du = 2.0 * arc.mConic[0] * u + arc.mConic[1] * v + arc.mConic[3];
            const double dv = arc.mConic[1] * u + 2.0 * arc.mConic[2] * v + arc.mConic[4];
            return du * eu + dv * ev;
        }

        static double implicitSecondDerivative(const ShadowDiscontinuityArc& arc, double eu, double ev)
        {
            if (arc.mType == ShadowDiscontinuityArc::SEGMENT)
            {
                return 0.0;
            }

            return arc.mConic[0] * eu * eu + arc.mConic[1] * eu * ev + arc.mConic[2] * ev * ev;
        }

        // an arc is a duplicate when it runs the same facets between the same endpoints, either way round
        static bool isDuplicate(const ShadowDiscontinuityArc& candidate, const std::vector<ShadowDiscontinuityArc>& output, double eps)
        {
            for (const ShadowDiscontinuityArc& other : output)
            {
                if (other.mFacetCount != candidate.mFacetCount)
                {
                    continue;
                }

                bool sameFacets = true;
                for (size_t i = 0; i < candidate.mFacetCount && sameFacets; i++)
                {
                    bool found = false;
                    for (size_t o = 0; o < other.mFacetCount && !found; o++)
                    {
                        found = other.mFacets[o] == candidate.mFacets[i];
                    }
                    sameFacets = found;
                }

                if (!sameFacets)
                {
                    continue;
                }

                const double d00 = (other.mP0 - candidate.mP0).getNorm();
                const double d11 = (other.mP1 - candidate.mP1).getNorm();
                const double d01 = (other.mP0 - candidate.mP1).getNorm();
                const double d10 = (other.mP1 - candidate.mP0).getNorm();

                if ((d00 <= eps && d11 <= eps) || (d01 <= eps && d10 <= eps))
                {
                    return true;
                }
            }

            return false;
        }

        /** @brief Split a raw arc where it crosses the receiver polygon boundary, keep the pieces running inside, append them.

        The crossing parameters are read off the implicit curve against the polygon edge lines (a
        quadratic per edge, a linear for a segment), turned back into values of the family parameter
        by eventParameterAt, and the ranges between them tested by their middle point. The pieces keep
        the implicit form and the parameterization of the raw arc: only the endpoints and the
        parameter range change. */
        static void clipAndEmit(const ShadowDiscontinuityArc& raw,
                                const std::vector<MathVector2d>& polygon,
                                const ReceiverFrame& frame,
                                std::vector<ShadowDiscontinuityArc>& output)
        {
            double tLow = std::min(raw.mT0, raw.mT1);
            double tHigh = std::max(raw.mT0, raw.mT1);
            const double span = tHigh - tLow;

            if (span <= 0.0 || !std::isfinite(span))
            {
                return;
            }

            const double tEps = span * 1e-6;
            std::vector<MathVector2d> params;
            params.push_back(MathVector2d(tLow, 0.0));
            // the range end as a sentinel, so the piece running from the last crossing to the end of the arc is tested too
            params.push_back(MathVector2d(tHigh, 0.0));

            const size_t count = polygon.size();
            for (size_t e = 0; e < count; e++)
            {
                const MathVector2d& v0 = polygon[e];
                const MathVector2d& v1 = polygon[(e + 1) % count];

                const double eu = v1.x - v0.x;
                const double ev = v1.y - v0.y;

                // points of the boundary line: p(lambda) = v0 + lambda (v1 - v0); f(p(lambda)) is quadratic
                std::vector<double> roots;
                solveQuadratic(implicitSecondDerivative(raw, eu, ev),
                               implicitDerivative(raw, v0.x, v0.y, eu, ev),
                               implicitAt(raw, v0.x, v0.y),
                               0.0, 1.0, roots);

                for (double lambda : roots)
                {
                    MathVector2d crossing(v0.x + lambda * eu, v0.y + lambda * ev);
                    double t = raw.eventParameterAt(crossing);

                    if (std::isnan(t))
                    {
                        continue;
                    }

                    if (t > tLow + tEps && t < tHigh - tEps)
                    {
                        params.push_back(MathVector2d(t, 0.0));
                    }
                }
            }

            std::sort(params.begin(), params.end(), [](const MathVector2d& a, const MathVector2d& b) { return a.x < b.x; });
            params.erase(std::unique(params.begin(), params.end(),
                                     [tEps](const MathVector2d& a, const MathVector2d& b) { return std::fabs(a.x - b.x) <= tEps; }),
                         params.end());

            for (size_t i = 0; i + 1 < params.size(); i++)
            {
                const double ta = params[i].x;
                const double tb = params[i + 1].x;

                MathVector2d middle = raw.pointAt(0.5 * (ta + tb));
                if (std::isnan(middle.x) || std::isnan(middle.y))
                {
                    continue;
                }

                if (!insidePolygon(polygon, frame, middle))
                {
                    continue;
                }

                MathVector2d pa = (ta == tLow) ? (raw.mT0 <= raw.mT1 ? raw.mP0 : raw.mP1) : raw.pointAt(ta);
                MathVector2d pb = (tb == tHigh) ? (raw.mT0 <= raw.mT1 ? raw.mP1 : raw.mP0) : raw.pointAt(tb);

                if (std::isnan(pa.x) || std::isnan(pa.y) || std::isnan(pb.x) || std::isnan(pb.y))
                {
                    continue;
                }

                MathVector2d chord = pb - pa;
                if (chord.getNorm() <= 10.0 * frame.mEps)
                {
                    continue;
                }

                ShadowDiscontinuityArc piece = raw;
                piece.mP0 = pa;
                piece.mP1 = pb;
                piece.mT0 = ta;
                piece.mT1 = tb;

                if (piece.mType == ShadowDiscontinuityArc::SEGMENT)
                {
                    // the endpoints moved: refit the line through them
                    MathVector3d line = homogeneousCross(pa, pb);
                    piece.mLine[0] = line.x;
                    piece.mLine[1] = line.y;
                    piece.mLine[2] = line.z;
                }

                if (isDuplicate(piece, output, 10.0 * frame.mEps))
                {
                    continue;
                }


                output.push_back(piece);
            }
        }

        // the implicit plane-coordinate line through the two endpoints of a segment arc
        static void fillSegmentLine(ShadowDiscontinuityArc& arc)
        {
            MathVector3d line = homogeneousCross(arc.mP0, arc.mP1);
            arc.mLine[0] = line.x;
            arc.mLine[1] = line.y;
            arc.mLine[2] = line.z;
        }

        /** @brief Whether two lines meet at a finite point, read on the Klein form with a relative tolerance.

        The Klein form vanishes for coplanar pairs, and parallel lines are coplanar: they meet at
        infinity, not at an event vertex, so the direction cross product must not vanish as well.
        Coincident lines have no direction cross either, and are left to the caller. */
        static bool linesMeet(const MathPlucker6d& a, const MathPlucker6d& b)
        {
            const double meetEpsilon = 1e-7;

            MathVector3d da = a.getDirection();
            MathVector3d db = b.getDirection();

            MathVector3d cross = MathVector3d::cross(da, db);
            double d2 = da.dot(da) * db.dot(db);
            if (d2 <= 0.0 || cross.dot(cross) <= 1e-10 * d2)
            {
                return false;
            }

            double k = std::fabs(a.dot(b));
            double n = a.scalarDot(a) * b.scalarDot(b);

            if (n <= 0.0)
            {
                return false;
            }

            return k <= meetEpsilon * MathArithmetic<double>::getSqrt(n);
        }
    };

    template<class P, class S>
    void ShadowDiscontinuityTracer<P, S>::traceLeafPolytope(PluckerPolyhedron<P>* aPolyhedron,
                                                            PluckerPolytope<P>* aPolytope,
                                                            const VisibilityEventRegistry& aRegistry,
                                                            const GeometryConvexPolygon& aSource,
                                                            const GeometryConvexPolygon& aReceiver,
                                                            S aTolerance,
                                                            std::vector<ShadowDiscontinuityArc>& anOutput)
    {
        const MathPlane3d sourcePlane = aSource.getPlane();
        const MathPlane3d receiverPlane = aReceiver.getPlane();

        const ReceiverFrame frame = makeFrame(aReceiver);
        if (frame.mEps <= 0.0)
        {
            return;
        }

        const std::vector<MathVector2d> polygon = projectPolygon(frame, aReceiver);

        const std::set<std::pair<size_t, size_t> >& edges = aPolytope->getEdges();
        const std::set<std::pair<size_t, size_t> >& realEdges = aPolytope->getEdgesIntersectingQuadric();


        for (const std::pair<size_t, size_t>& edge : edges)
        {
            if constexpr (!MathPluckerCellTraits<P>::cAllLinesReal)
            {
                // only the edges crossing the quadric hold real lines, and real lines are what traces a swath on the receiver
                if (realEdges.find(edge) == realEdges.end())
                {
                    continue;
                }
            }

            std::vector<size_t> common;
            MathCombinatorial::getCommonFacets(aPolyhedron->getFacetsDescription(edge.first), aPolyhedron->getFacetsDescription(edge.second), common);

            if (common.size() != MathPluckerCellTraits<P>::cCommonFacetsForEdge)
            {
                // more shared facets: the two vertices stand at the same line, a zero length family
                continue;
            }

            /* The two lines bounding the family on the receiver: the crossings of the edge with the
            Plucker quadric, or, on a family whose every line is real already, the corner rays of
            the cell itself. */
            std::vector<P> endpoints;

            if constexpr (MathPluckerCellTraits<P>::cAllLinesReal)
            {
                endpoints.push_back(aPolyhedron->get(edge.first));
                endpoints.push_back(aPolyhedron->get(edge.second));
            }
            else
            {
                std::vector<P> crossings;
                if (!MathGeometry::findPluckerEdgeWithQuadricIntersection<P, S>(aPolyhedron->get(edge.first),
                                                                                aPolyhedron->get(edge.second),
                                                                                aPolyhedron->getQuadricRelativePosition(edge.first),
                                                                                aPolyhedron->getQuadricRelativePosition(edge.second),
                                                                                crossings, aTolerance))
                {
                    continue;
                }

                if (crossings.size() < 2)
                {
                    // a tangent edge: the arc has no length to trace
                    continue;
                }

                endpoints.swap(crossings);
            }

            MathVector3d rcv0 = MathGeometry::getBackTo3D<P>(endpoints[0], sourcePlane, receiverPlane).second;
            MathVector3d rcv1 = MathGeometry::getBackTo3D<P>(endpoints[1], sourcePlane, receiverPlane).second;

            if (!std::isfinite(rcv0.x) || !std::isfinite(rcv0.y) || !std::isfinite(rcv0.z)
                || !std::isfinite(rcv1.x) || !std::isfinite(rcv1.y) || !std::isfinite(rcv1.z))
            {
                continue;
            }

            // a stabbing line grazing the receiver plane pierces it anywhere: a trace this far out is an artifact
            MathVector3d reach0 = rcv0 - frame.mOrigin;
            MathVector3d reach1 = rcv1 - frame.mOrigin;
            if (reach0.getNorm() > 100.0 * frame.mSize || reach1.getNorm() > 100.0 * frame.mSize)
            {
                continue;
            }

            MathVector2d p0 = project(frame, rcv0);
            MathVector2d p1 = project(frame, rcv1);

            MathVector2d chord = p1 - p0;
            if (chord.getNorm() <= 10.0 * frame.mEps)
            {
                continue;
            }

            ShadowDiscontinuityArc arc;
            arc.mType = ShadowDiscontinuityArc::SEGMENT;
            arc.mP0 = p0;
            arc.mP1 = p1;
            arc.mFacetCount = common.size();
            for (size_t f = 0; f < common.size(); f++)
            {
                arc.mFacets[f] = common[f];
            }
            arc.mHasOccluderFacet = false;
            arc.mOrigin = frame.mOrigin;
            arc.mE1 = frame.mE1;
            arc.mE2 = frame.mE2;
            arc.mT0 = 0.0;
            arc.mT1 = 1.0;

            if constexpr (MathPluckerCellTraits<P>::cAllLinesReal)
            {
                /* Directional light: the edge runs on the shadow line of one carrier edge, and the
                swath is the plane of that edge and the direction of the light. The traces of the
                source boundaries on the receiver are not shadow discontinuities, and only the
                silhouette edges cast shadows. */
                const HyperplaneEventInfo* info = aRegistry.find(common[0]);
                if (info == nullptr || info->mTouchesSynthetic || info->mKind != HyperplaneEventInfo::SilhouetteEdge)
                {
                    continue;
                }

                arc.mHasOccluderFacet = true;
                fillSegmentLine(arc);
                clipAndEmit(arc, polygon, frame, anOutput);
                continue;
            }

            // the three event lines the family lies on
            MathPlucker6d L[3];
            MathVector3d x[3];        // anchors
            MathVector3d u[3];        // normalized directions
            MathVector3d dStored[3];  // directions and moments as stored: the scales a row of F may bring along cancel out
            MathVector3d mStored[3];
            const HyperplaneEventInfo* infos[3];

            for (size_t f = 0; f < 3; f++)
            {
                L[f] = convert<MathPlucker6d>(aPolyhedron->get(common[f]));
                dStored[f] = L[f].getDirection();
                mStored[f] = L[f].getLocation();

                x[f] = MathGeometry::getAnchorPoint<MathPlucker6d>(L[f]);
                u[f] = dStored[f];
                u[f].normalize();

                infos[f] = aRegistry.find(common[f]);

                // an unregistered line may be anything but a registered source edge: it counts as an occluder
                if (infos[f] == nullptr || infos[f]->mKind == HyperplaneEventInfo::SilhouetteEdge)
                {
                    arc.mHasOccluderFacet = true;
                }
            }

            const bool meet[3] = { linesMeet(L[0], L[1]), linesMeet(L[0], L[2]), linesMeet(L[1], L[2]) };
            const int meetCount = (meet[0] ? 1 : 0) + (meet[1] ? 1 : 0) + (meet[2] ? 1 : 0);


            if (meetCount >= 2)
            {
                /* Lines meeting two by two or more bound no swath worth tracing here: the family
                collapses (three coplanar lines) or turns degenerate; the event is left out. */
                continue;
            }

            if (meetCount == 1)
            {
                // a VE event: the swath is the plane of the meeting vertex and the remaining line
                const int i = meet[2] ? 1 : 0;
                const int j = meet[0] ? 1 : 2;
                const int k = 3 - i - j;

                MathVector3d vertex;
                if (!intersectPoint(x[i], u[i], x[j], u[j], vertex))
                {
                    // nearly parallel: the meeting is read too ill conditioned to trust
                    continue;
                }

                // the vertex must stand on both lines: a grazing reading sends the swath plane astray
                double vScale = std::max(frame.mSize, (vertex - frame.mOrigin).getNorm());
                if (MathVector3d::cross(u[i], vertex - x[i]).getNorm() > 1e-5 * vScale
                    || MathVector3d::cross(u[j], vertex - x[j]).getNorm() > 1e-5 * vScale)
                {
                    continue;
                }

                if (infos[i] != nullptr && infos[j] != nullptr
                    && infos[i]->mKind == HyperplaneEventInfo::SourceEdge && infos[j]->mKind == HyperplaneEventInfo::SourceEdge
                    && infos[i]->mPolygonIndex == infos[j]->mPolygonIndex
                    && (infos[i]->mTouchesSynthetic || infos[j]->mTouchesSynthetic))
                {
                    // a vertex of the padding that makes a degenerated source a polygon: an artifact, not an occluder corner
                    continue;
                }

                MathVector3d swathNormal = MathVector3d::cross(u[k], vertex - x[k]);
                double sn = swathNormal.normalize();
                if (sn <= 1e-9)
                {
                    // the vertex runs on the remaining line: the family is a whole pencil, not a swath
                    continue;
                }

                double inclination = swathNormal.dot(frame.mNormal);
                if (1.0 - inclination * inclination <= 1e-12)
                {
                    // a swath parallel to the receiver: the trace runs to infinity
                    continue;
                }

                arc.mType = ShadowDiscontinuityArc::SEGMENT;
                arc.mT0 = 0.0;
                arc.mT1 = 1.0;
                fillSegmentLine(arc);
            }
            else
            {
                // an EEE event: the transversals sweep a regulus, its conic on the receiver plane
                double coeffs[6];
                if (!buildConic(dStored, mStored, frame, coeffs))
                {
                    continue;
                }

                arc.mType = ShadowDiscontinuityArc::CONIC;
                for (int c = 0; c < 6; c++)
                {
                    arc.mConic[c] = coeffs[c];
                }

                // the endpoints of the arc must stand on its own curve, else the reading of the event is unreliable
                if (std::fabs(implicitAt(arc, arc.mP0.x, arc.mP0.y)) > 1e-2 || std::fabs(implicitAt(arc, arc.mP1.x, arc.mP1.y)) > 1e-2)
                {
                    continue;
                }

                // parameterize the family on whichever event line reads both endpoint parameters
                MathPlucker6d eL[2] = { convert<MathPlucker6d>(endpoints[0]), convert<MathPlucker6d>(endpoints[1]) };
                MathVector3d ex[2];
                MathVector3d eu[2];
                for (int q = 0; q < 2; q++)
                {
                    ex[q] = MathGeometry::getAnchorPoint<MathPlucker6d>(eL[q]);
                    eu[q] = eL[q].getDirection();
                    eu[q].normalize();
                }

                bool parametrized = false;
                for (int par = 0; par < 3 && !parametrized; par++)
                {
                    double t0 = 0.0;
                    double t1 = 0.0;
                    if (!ShadowDiscontinuityArc::meetOnLine(x[par], u[par], ex[0], eu[0], t0)
                        || !ShadowDiscontinuityArc::meetOnLine(x[par], u[par], ex[1], eu[1], t1)
                        || t0 == t1)
                    {
                        continue;
                    }

                    int others[2];
                    for (int o = 0, q = 0; q < 3; q++)
                    {
                        if (q != par)
                        {
                            others[o++] = q;
                        }
                    }

                    arc.mA[0] = x[par];
                    arc.mU[0] = u[par];
                    arc.mA[1] = x[others[0]];
                    arc.mU[1] = u[others[0]];
                    arc.mA[2] = x[others[1]];
                    arc.mU[2] = u[others[1]];
                    arc.mT0 = t0;
                    arc.mT1 = t1;
                    parametrized = true;
                }

                if (!parametrized)
                {
                    continue;
                }
            }

            // bounds of the sources themselves trace no shadow: an arc needs a silhouette edge among its facets
            if (!arc.mHasOccluderFacet)
            {
                continue;
            }

            clipAndEmit(arc, polygon, frame, anOutput);
        }

    }
}
