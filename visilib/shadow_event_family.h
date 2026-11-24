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
#include <utility>
#include <vector>

#include "math_plucker_6.h"
#include "math_polynomial.h"
#include "math_vector_3.h"

namespace visilib
{
    /** @brief A half space of a cell, read in double precision: the cell lies where mSign * (L . h) <= 0.

    mSign is 1 when the cell lies where the Klein form with the hyperplane is negative, -1 where it is positive, and 0
    when the whole cell lies in the hyperplane (an equality), as PluckerPolytope::CellFacet gives it. The hyperplane is
    normalized by its six dimensional norm.*/
    struct ShadowCellFacet
    {
        MathVector3d mDirection;
        MathVector3d mLocation;
        double mSign;
        size_t mGroup;  /**< @brief The collinear group of the edge the hyperplane carries (see VisibilityEventRegistry), or -1*/
    };

    /** @brief A one parameter family of lines, the critical lines of a visibility event, and the part of it a cell holds.

    The lines are given as polynomials of the parameter t, in Plucker coordinates up to a common factor:

        L(t) = ( sum_k mDirection[k] t^k , sum_k mLocation[k] t^k )

    The Klein form of L(t) with any fixed line is then a polynomial of t, whose roots are where L(t) crosses the
    hyperplane of that line. The three families of the shadow events are all of second degree at most:

    - the vertex pencil, the lines through a vertex v meeting the segment [A, B], through p(t) = A + t (B - A):
      L(t) = (p(t) - v, v x p(t)), of first degree.
    - the directional family, the rays of direction d meeting [A, B]: L(t) = (d, p(t) x d), of first degree.
    - the transversals of three lines, through the point a(t) = A1 + t E1 of the first one. The plane of a(t) and of
      the second line (P2, E2) has the normal n(t) = E2 x (a(t) - P2), of first degree, and meets the third line
      (P3, E3) at c(t) = P3 + E3 N(t) / D(t), with N = n . (P2 - P3) and D = n . E3, both of first degree. With
      C = P3 D + E3 N, the transversal is L(t) = (C - a D, a x C) / D: the numerator is of second degree, and the
      common factor 1 / D is dropped (its sign is taken back by the orientation, see getLine).

    The lines a family holds as real visibility events are those oriented from the light to the receiver: getLine
    orients L(t) along mOrientation, the vector from the light to the receiver.*/
    class ShadowEventFamily
    {
    public:
        /** @brief The lines through aVertex meeting the segment [aBegin, anEnd], for t in [0, 1]. False when the vertex is on the line of the segment.*/
        bool buildVertexPencil(const MathVector3d& aVertex, const MathVector3d& aBegin, const MathVector3d& anEnd, const MathVector3d& anOrientation)
        {
            clear(anOrientation);
            const MathVector3d e = anEnd - aBegin;
            if (MathVector3d::cross(e, aBegin - aVertex).getSquaredNorm() <= cDegenerate * e.getSquaredNorm() * (aBegin - aVertex).getSquaredNorm())
            {
                return false;
            }
            mDirection = { aBegin - aVertex, e };
            mLocation = { MathVector3d::cross(aVertex, aBegin), MathVector3d::cross(aVertex, e) };
            return true;
        }

        /** @brief The rays of direction aDirection meeting the segment [aBegin, anEnd], for t in [0, 1]. False when the segment is along the rays.*/
        bool buildDirectional(const MathVector3d& aDirection, const MathVector3d& aBegin, const MathVector3d& anEnd)
        {
            clear(aDirection);
            const MathVector3d e = anEnd - aBegin;
            if (MathVector3d::cross(e, aDirection).getSquaredNorm() <= cDegenerate * e.getSquaredNorm() * aDirection.getSquaredNorm())
            {
                return false;
            }
            mDirection = { aDirection, MathVector3d::Zero() };
            mLocation = { MathVector3d::cross(aBegin, aDirection), MathVector3d::cross(e, aDirection) };
            return true;
        }

        /** @brief The transversals of the three segments [b[i], e[i]], through the points of the first one, for t in [0, 1].

        The second and the third segments bound the family too: the lines must meet them between their endpoints, which
        the breakpoints of their meeting parameters and the support test of isFree enforce. False when the parameterization
        degenerates (the third line parallel to the plane of the second one and of every point of the first).*/
        bool buildTransversals(const MathVector3d b[3], const MathVector3d e[3], const MathVector3d& anOrientation)
        {
            clear(anOrientation);
            const MathVector3d A1 = b[0];
            const MathVector3d E1 = e[0] - b[0];
            const MathVector3d P2 = b[1];
            const MathVector3d E2 = e[1] - b[1];
            const MathVector3d P3 = b[2];
            const MathVector3d E3 = e[2] - b[2];

            // n(t) = n0 + t n1, N(t) = N0 + t N1, D(t) = D0 + t D1
            const MathVector3d n0 = MathVector3d::cross(E2, A1 - P2);
            const MathVector3d n1 = MathVector3d::cross(E2, E1);
            const double N0 = n0.dot(P2 - P3), N1 = n1.dot(P2 - P3);
            const double D0 = n0.dot(E3), D1 = n1.dot(E3);

            const double scale = E1.getNorm() * E2.getNorm() * E3.getNorm();
            if (std::max(std::fabs(D0), std::fabs(D1)) <= cDegenerate * scale)
            {
                return false;
            }

            // C(t) = P3 D(t) + E3 N(t) = C0 + t C1 ; a(t) = A1 + t E1
            const MathVector3d C0 = P3 * D0 + E3 * N0;
            const MathVector3d C1 = P3 * D1 + E3 * N1;

            // direction C - a D, location a x C, both of second degree
            mDirection = { C0 - A1 * D0, C1 - A1 * D1 - E1 * D0, E1 * (-D1) };
            mLocation = { MathVector3d::cross(A1, C0), MathVector3d::cross(A1, C1) + MathVector3d::cross(E1, C0), MathVector3d::cross(E1, C1) };

            // the third segment is met at the parameter N / D: it changes of side of [0, 1] at the roots of N, N - D and D
            mBreakpoints.push_back(MathPolynomiald({ N0, N1 }));
            mBreakpoints.push_back(MathPolynomiald({ N0 - D0, N1 - D1 }));
            mBreakpoints.push_back(MathPolynomiald({ D0, D1 }));

            // the second segment, by the plane of a(t) and of the third line: n'(t) = E3 x (a(t) - P3)
            const MathVector3d m0 = MathVector3d::cross(E3, A1 - P3);
            const MathVector3d m1 = MathVector3d::cross(E3, E1);
            const double M0 = m0.dot(P3 - P2), M1 = m1.dot(P3 - P2);
            const double F0 = m0.dot(E2), F1 = m1.dot(E2);
            mBreakpoints.push_back(MathPolynomiald({ M0, M1 }));
            mBreakpoints.push_back(MathPolynomiald({ M0 - F0, M1 - F1 }));
            mBreakpoints.push_back(MathPolynomiald({ F0, F1 }));

            mSupports.push_back(std::make_pair(b[1], e[1]));
            mSupports.push_back(std::make_pair(b[2], e[2]));
            return true;
        }

        /** @brief The line L(t), oriented from the light to the receiver. False when it is orthogonal to that direction (or null).*/
        bool getLine(double t, MathPlucker6d& aLine) const
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
            const double orientation = d.dot(mOrientation);
            const double dn = d.getNorm();
            if (!(dn > 0.0) || std::fabs(orientation) <= cDegenerate * dn * mOrientation.getNorm())
            {
                return false;
            }
            if (orientation < 0.0)
            {
                d *= -1.0;
                m *= -1.0;
            }
            // unit direction: the Klein forms read as distances times sines
            d *= 1.0 / dn;
            m *= 1.0 / dn;
            aLine = MathPlucker6d(d.x, d.y, d.z, m.x, m.y, m.z);
            return true;
        }

        /** @brief The Klein form of L(t) with the line (aDirection, aLocation), as a polynomial of t.*/
        MathPolynomiald getKleinPolynomial(const MathVector3d& aDirection, const MathVector3d& aLocation) const
        {
            std::vector<double> c(mDirection.size());
            for (size_t k = 0; k < mDirection.size(); k++)
            {
                c[k] = mDirection[k].dot(aLocation) + mLocation[k].dot(aDirection);
            }
            return MathPolynomiald(c);
        }

        /** @brief Whether the line L(t) lies in the cell described by someFacets and meets the support segments.

        aTolerance is a distance: a line may stand that far outside a facet, or beyond the end of a segment (relative to its
        length), and still be taken as in.*/
        bool isFree(double t, const std::vector<ShadowCellFacet>& someFacets, double aTolerance) const
        {
            MathPlucker6d line;
            if (!getLine(t, line))
            {
                return false;
            }
            const double norm = MathArithmetic<double>::getSqrt(line.getDirection().getSquaredNorm() + line.getLocation().getSquaredNorm());
            for (const ShadowCellFacet& facet : someFacets)
            {
                const double value = (line.getDirection().dot(facet.mLocation) + line.getLocation().dot(facet.mDirection)) / norm;
                if (facet.mSign == 0.0)
                {
                    if (std::fabs(value) > aTolerance)
                    {
                        return false;
                    }
                }
                else if (facet.mSign * value > aTolerance)
                {
                    return false;
                }
            }
            for (const auto& segment : mSupports)
            {
                const MathPlucker6d support(segment.first, segment.second);
                MathVector3d meeting;
                if (!support.getMeetingPoint(line, meeting, 1e-24))
                {
                    return false;
                }
                const MathVector3d e = segment.second - segment.first;
                const double s = (meeting - segment.first).dot(e) / e.getSquaredNorm();
                if (s < -aTolerance || s > 1.0 + aTolerance)
                {
                    return false;
                }
            }
            return true;
        }

        /** @brief The ranges of t in [0, 1] where L(t) is free in the cell: append them to intervals, sorted and disjoint.

        The breakpoints are the roots of the Klein polynomials of the facets (where L(t) crosses a facet), of the
        orientation (where L(t) turns away from the receiver) and of the support parameters. Between two consecutive
        breakpoints nothing changes, and the middle of the range decides for all of it (see isFree). Consecutive free
        ranges are joined.*/
        void getFreeIntervals(const std::vector<ShadowCellFacet>& someFacets, double aTolerance, std::vector<std::pair<double, double>>& intervals) const
        {
            const double relative = 1e-12;
            std::vector<double> breakpoints = { 0.0, 1.0 };
            for (const ShadowCellFacet& facet : someFacets)
            {
                getKleinPolynomial(facet.mDirection, facet.mLocation).findRoots(0.0, 1.0, relative, breakpoints);
            }
            std::vector<double> orientation(mDirection.size());
            for (size_t k = 0; k < mDirection.size(); k++)
            {
                orientation[k] = mDirection[k].dot(mOrientation);
            }
            MathPolynomiald(orientation).findRoots(0.0, 1.0, relative, breakpoints);
            for (const MathPolynomiald& p : mBreakpoints)
            {
                p.findRoots(0.0, 1.0, relative, breakpoints);
            }
            std::sort(breakpoints.begin(), breakpoints.end());

            const double minimumLength = 1e-9;
            std::vector<std::pair<double, double>> free;
            for (size_t i = 0; i + 1 < breakpoints.size(); i++)
            {
                const double t0 = breakpoints[i];
                const double t1 = breakpoints[i + 1];
                if (t1 - t0 < minimumLength)
                {
                    continue;
                }
                if (!isFree(0.5 * (t0 + t1), someFacets, aTolerance))
                {
                    continue;
                }
                if (!free.empty() && t0 - free.back().second < minimumLength)
                {
                    free.back().second = t1;
                }
                else
                {
                    free.push_back(std::make_pair(t0, t1));
                }
            }
            intervals.insert(intervals.end(), free.begin(), free.end());
        }

        const std::vector<MathVector3d>& getDirectionCoefficients() const { return mDirection; }
        const std::vector<MathVector3d>& getLocationCoefficients() const { return mLocation; }

    private:
        void clear(const MathVector3d& anOrientation)
        {
            mDirection.clear();
            mLocation.clear();
            mBreakpoints.clear();
            mSupports.clear();
            mOrientation = anOrientation;
        }

        static constexpr double cDegenerate = 1e-18; /**< @brief A squared sine below which two directions are taken as parallel*/

        std::vector<MathVector3d> mDirection;                            /**< @brief The coefficients of the direction of L(t)*/
        std::vector<MathVector3d> mLocation;                             /**< @brief The coefficients of the location of L(t)*/
        std::vector<MathPolynomiald> mBreakpoints;                       /**< @brief Polynomials whose roots bound where the supports are met*/
        std::vector<std::pair<MathVector3d, MathVector3d>> mSupports;    /**< @brief The segments L(t) must meet between their endpoints*/
        MathVector3d mOrientation;                                       /**< @brief The direction from the light to the receiver*/
    };
}
