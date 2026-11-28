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
#include <cstdint>
#include <type_traits>
#include <vector>

#include "geometry_position_type.h"
#include "math_plucker_6.h"
#include "visilib_core.h"

namespace visilib
{
    /** @brief Whether the arithmetic S is exact: its sums, products and quotients carry no rounding, and its signs are
    true. Only then does Simulation of Simplicity make sense: a perturbation can only break ties that are real zeros.*/
    template<class S>
    struct MathIsExact
    {
        static constexpr bool value = false;
    };

#ifdef ENABLE_GMP
    template<>
    struct MathIsExact<MathGmpRational>
    {
        static constexpr bool value = true;
    };
#endif
#if defined(ENABLE_GMP) && defined(ENABLE_MPFI)
    /* The lazy rational has the values of GMP rational: its comparisons are exact, only decided by intervals first.*/
    template<>
    struct MathIsExact<MathLazyRational>
    {
        static constexpr bool value = true;
    };
#endif

    /** @brief Whether the double sign filter (MathSignFilter) applies to the arithmetic S: its conversion to double
    must round the exact value once. The lazy rational converts to the middle of an interval, which can be wide, and
    filters its comparisons itself.*/
    template<class S>
    struct MathHasSignFilter
    {
        static constexpr bool value = false;
    };
#ifdef ENABLE_GMP
    template<>
    struct MathHasSignFilter<MathGmpRational>
    {
        static constexpr bool value = true;
    };
#endif

    /** @brief Whether Simulation of Simplicity applies to the Plucker points P of the arithmetic S: an exact arithmetic,
    and the general lines (the directional family re-projects its moment after each sum, which a polynomial of epsilon
    cannot carry through).*/
    template<class P, class S>
    struct MathSoSSupported
    {
        static constexpr bool value = MathIsExact<S>::value && std::is_same<P, MathPlucker6<S>>::value;
    };

    /** @brief A floating-point filter for the sign of the Plucker product of two exact points.

    The exact coordinates of the deep vertices of the CSG grow to thousands of bits, and an exact product of them costs
    far more than its sign is worth: almost every sign is far from zero. The filter computes the product of the double
    approximations with a bound of its error, and tells the sign only when the bound cannot change it; otherwise (0)
    the exact product decides. The bound:
    - each coordinate is rounded once (mpq_get_d truncates: a relative error below 2^-52), or flushed to zero below
      the normal range (an absolute error below 2^-1021);
    - the six products and their sum in double add at most 7 roundings of 2^-53, relative to the sum of the absolute
      products.
    So |exact - computed| <= 16 * 2^-52 * sum |a_i| |b_j| + 2^-1000 * (sum |a_i| + sum |b_j|) * (1 + max), well inside
    the bound below. Infinite or not-a-number approximations (a coordinate beyond the double range) decide nothing.*/
    struct MathSignFilter
    {
        /** @brief The approximation of a point: its coordinates rounded to double.*/
        template<class S>
        static MathPlucker6d approximate(const MathPlucker6<S>& aPoint)
        {
            return convert<MathPlucker6d>(aPoint);
        }

        /** @brief The sign of the Plucker product of the exact points that a and b approximate, or 0 when the
        approximations cannot tell it. With aThreshold, the sign is told only when the absolute value of the product is
        also above it (a tolerance that calls the smaller values zero).*/
        static int sign(const MathPlucker6d& a, const MathPlucker6d& b, double aThreshold = 0.0)
        {
            const double value = a.dot(b);
            const MathVector3d& ad = a.getDirection();
            const MathVector3d& al = a.getLocation();
            const MathVector3d& bd = b.getDirection();
            const MathVector3d& bl = b.getLocation();
            const double magnitude = std::fabs(ad.x * bl.x) + std::fabs(ad.y * bl.y) + std::fabs(ad.z * bl.z)
                                   + std::fabs(al.x * bd.x) + std::fabs(al.y * bd.y) + std::fabs(al.z * bd.z);
            const double sizeA = std::fabs(ad.x) + std::fabs(ad.y) + std::fabs(ad.z) + std::fabs(al.x) + std::fabs(al.y) + std::fabs(al.z);
            const double sizeB = std::fabs(bd.x) + std::fabs(bd.y) + std::fabs(bd.z) + std::fabs(bl.x) + std::fabs(bl.y) + std::fabs(bl.z);
            // beyond the double range, or so small that the products underflow: the filter does not decide
            if (!std::isfinite(value) || !std::isfinite(magnitude) || !std::isfinite(sizeA) || !std::isfinite(sizeB) || magnitude < 1.0e-250)
            {
                return 0;
            }
            const double bound = 4.0e-15 * magnitude + 1.0e-290 * (sizeA + sizeB) * (1.0 + sizeA + sizeB);
            const double limit = bound + aThreshold;
            return value > limit ? 1 : (value < -limit ? -1 : 0);
        }
    };

    /** @brief The exact sign of a value: -1, 0 or +1*/
    template<class S>
    inline int getExactSign(const S& x)
    {
        return (x > S(0)) ? 1 : ((x < S(0)) ? -1 : 0);
    }

    /** @brief A Plucker point of the symbolically perturbed scene: a polynomial of the infinitesimal epsilon,

        v(epsilon) = c[0] + c[1] epsilon + c[2] epsilon^2 + ...

    c[0] is the point the polyhedron stores; the higher terms say where the perturbation took it. A point of the input
    (a vertex built from the source and the receiver) has no higher term. The leading term c[0] is never zero: a
    polynomial whose first terms vanish is divided by the power of epsilon they make, a positive factor, which a
    projective point does not see.*/
    template<class P>
    struct MathSymbolicPoint
    {
        std::vector<P> c;

        MathSymbolicPoint() {}
        explicit MathSymbolicPoint(const P& aLeading) : c(1, aLeading) {}

        const P& leading() const { return c[0]; }
        size_t size() const { return c.size(); }

        /** @brief Multiply every coefficient by a positive factor: the same projective point*/
        template<class S>
        void scale(const S& aFactor)
        {
            for (P& p : c)
            {
                p *= aFactor;
            }
        }
    };

    /** @brief Simulation of Simplicity for the Plucker CSG (Edelsbrunner and Mucke), in an exact arithmetic.

    The silhouette hyperplanes are perturbed: the hyperplane of index i becomes h_i + epsilon u_i, u_i a small integer
    vector derived from i, different from one hyperplane to the next. The vertices the splits build from them carry the
    perturbation as the higher terms of their symbolic point. Every decision is then the sign of a polynomial of
    epsilon, read from its first coefficient that does not vanish: the decisions are the ones of a single scene in
    general position, epsilon being as small as needed, and the CSG only meets the generic case.

    A polynomial that vanishes to every order kept is not an accident the perturbation can break: the vertex lies on
    the hyperplane (or on the quadric) by construction, as every combination of two lines through one vertex of the
    light is again such a line. It is answered as an incidence, as the facet descriptions answer the recorded ones.*/
    template<class P, class S>
    class MathSoS
    {
    public:
        /** @brief The number of terms of epsilon kept beyond the leading one.

        Each term kept multiplies the cost of the crossings (a product of polynomials). On the test scenes every tie was
        broken at order 1, and none deeper; the second order is kept as a margin for the nested degeneracies, a tie on
        a point built from a broken tie. A polynomial that vanishes to that order is answered as an incidence.*/
        static constexpr size_t cOrder = 2;

        static bool isZero(const P& p)
        {
            return getExactSign(p.getDirection().x) == 0 && getExactSign(p.getDirection().y) == 0 && getExactSign(p.getDirection().z) == 0
                && getExactSign(p.getLocation().x) == 0 && getExactSign(p.getLocation().y) == 0 && getExactSign(p.getLocation().z) == 0;
        }

        /** @brief The perturbation direction of the hyperplane anIndex: six integers in [-7, 7], not all zero, a
        deterministic hash of the index.*/
        static P perturbationOf(size_t anIndex)
        {
            uint64_t x = uint64_t(anIndex) * 0x9E3779B97F4A7C15ull + 0x632BE59BD9B4E019ull;
            int c[6];
            bool allZero = true;
            for (int k = 0; k < 6; k++)
            {
                x ^= x >> 31;
                x *= 0xBF58476D1CE4E5B9ull;
                x ^= x >> 29;
                c[k] = int(x % 15) - 7;
                allZero &= c[k] == 0;
            }
            if (allZero)
            {
                c[0] = 1;
            }
            return P(S(c[0]), S(c[1]), S(c[2]), S(c[3]), S(c[4]), S(c[5]));
        }

        /** @brief The sign of the Klein product of the perturbed hyperplane h + epsilon u with the symbolic point v.

        The coefficient of epsilon^k is h . v_k + u . v_(k-1). anOrder returns the k that decided, or cOrder + 1 when
        every coefficient vanished: the return value is then 0, an incidence by construction.
        @param aPerturbation: u, null for a hyperplane that is not perturbed*/
        static int sign(const P& aPlane, const P* aPerturbation, const MathSymbolicPoint<P>& aPoint, size_t& anOrder)
        {
            for (size_t k = 0; k <= cOrder; k++)
            {
                S value(0);
                if (k < aPoint.size())
                {
                    value = aPlane.dot(aPoint.c[k]);
                }
                if (aPerturbation != nullptr && k >= 1 && k - 1 < aPoint.size())
                {
                    value = value + aPerturbation->dot(aPoint.c[k - 1]);
                }
                const int s = getExactSign(value);
                if (s != 0)
                {
                    anOrder = k;
                    return s;
                }
            }
            anOrder = cOrder + 1;
            return 0;
        }

        /** @brief The sign of the Plucker quadric at the symbolic point: the coefficient of epsilon^k of v . v is the sum
        of the v_a . v_b with a + b = k. 0 when every coefficient vanishes: the point is a real line by construction.*/
        static int quadricSign(const MathSymbolicPoint<P>& aPoint, size_t& anOrder)
        {
            for (size_t k = 0; k <= cOrder; k++)
            {
                S value(0);
                for (size_t a = 0; a <= k && a < aPoint.size(); a++)
                {
                    const size_t b = k - a;
                    if (b < aPoint.size())
                    {
                        value = value + aPoint.c[a].dot(aPoint.c[b]);
                    }
                }
                const int s = getExactSign(value);
                if (s != 0)
                {
                    anOrder = k;
                    return s;
                }
            }
            anOrder = cOrder + 1;
            return 0;
        }

        /** @brief The point where the edge from aNegative to aPositive crosses the perturbed hyperplane h + epsilon u.

        With s_n and s_p the polynomials of the hyperplane at the two ends (negative and positive), the point is
        s_p v_n - s_n v_p: both weights are positive, the point is on the segment between the ends, and it is division
        free. The product of polynomials is truncated at cOrder, and the leading zero coefficients, which a crossing at
        an end whose value vanished at order 0 leaves, are shifted out.*/
        static MathSymbolicPoint<P> crossing(const P& aPlane, const P* aPerturbation, const MathSymbolicPoint<P>& aNegative, const MathSymbolicPoint<P>& aPositive)
        {
            const std::vector<S> sn = values(aPlane, aPerturbation, aNegative);
            const std::vector<S> sp = values(aPlane, aPerturbation, aPositive);

            MathSymbolicPoint<P> result;
            for (size_t k = 0; k <= cOrder; k++)
            {
                P term = P::Zero();
                bool any = false;
                for (size_t a = 0; a <= k; a++)
                {
                    const size_t b = k - a;
                    if (b < aNegative.size() && getExactSign(sp[a]) != 0)
                    {
                        term += aNegative.c[b] * sp[a];
                        any = true;
                    }
                    if (b < aPositive.size() && getExactSign(sn[a]) != 0)
                    {
                        term -= aPositive.c[b] * sn[a];
                        any = true;
                    }
                }
                result.c.push_back(any ? term : P::Zero());
            }

            /* Bring the point back to the scale of its ends: divided by s_p - s_n at order 0, positive since the ends are on
            opposite sides, the leading term is the convex combination the interpolation of the rounded arithmetics gives,
            exactly. Without it the coordinates multiply at each split, and the parts of the solver that read the points
            in double (the representative line, the rays) lose them. The factor is positive: the same projective point. */
            const S weightSum = sp[0] - sn[0];
            if (getExactSign(weightSum) > 0)
            {
                const S factor = S(1) / weightSum;
                for (P& term : result.c)
                {
                    term *= factor;
                }
            }

            // the leading zeros: divide by the power of epsilon they make
            size_t first = 0;
            while (first < result.c.size() && isZero(result.c[first]))
            {
                first++;
            }
            if (first == result.c.size())
            {
                // the two ends are the same symbolic point: no crossing to build, keep the negative end
                return aNegative;
            }
            result.c.erase(result.c.begin(), result.c.begin() + first);
            // drop the trailing zeros, so that the points of the input stay with their single term
            while (result.c.size() > 1 && isZero(result.c.back()))
            {
                result.c.pop_back();
            }
            return result;
        }

        /** @brief Whether two symbolic points are the same, coefficient by coefficient*/
        static bool areEqual(const MathSymbolicPoint<P>& a, const MathSymbolicPoint<P>& b)
        {
            const size_t n = std::max(a.size(), b.size());
            for (size_t k = 0; k < n; k++)
            {
                const P pa = k < a.size() ? a.c[k] : P::Zero();
                const P pb = k < b.size() ? b.c[k] : P::Zero();
                if (!isZero(pa - pb))
                {
                    return false;
                }
            }
            return true;
        }

        /** @brief Whether the edge between two points of the polyhedron holds a real line, decided exactly.

        anEnd1 and anEnd2 are the positions of the ends relative to the quadric, as the polyhedron recorded them (an end
        on the quadric is a real line: the edge holds one). Ends on opposite sides: the edge crosses the quadric. Ends on
        the same side: the quadratic f(t) = a t^2 + 2 b t + c of the points v1 + t (v2 - v1) has two roots in the edge
        when its apex -b/a lies strictly inside (0, 1) and its value there, (a c - b^2) / a, has the opposite sign of
        the ends. A tangency, a discriminant exactly zero, is one real line and no interior: it counts as no crossing,
        and is reported through isTangent.*/
        static bool edgeCrossesQuadric(const P& v1, const P& v2, GeometryPositionType anEnd1, GeometryPositionType anEnd2, bool& isTangent)
        {
            isTangent = false;
            if (anEnd1 == ON_BOUNDARY || anEnd2 == ON_BOUNDARY)
            {
                return true;
            }
            if (anEnd1 != anEnd2)
            {
                return true;
            }

            const P f = v2 - v1;
            const S a = f.dot(f);
            const S b = f.dot(v1);
            const S c = v1.dot(v1);
            const int sa = getExactSign(a);
            if (sa == 0)
            {
                // f linear: monotone, and the ends are on the same side
                return false;
            }
            const S discriminant = b * b - a * c;
            const int sd = getExactSign(discriminant);
            if (sd < 0)
            {
                return false;
            }
            // the apex t* = -b / a strictly inside (0, 1): -b and a of the same sign, and |b| < |a|
            const S apexNumerator = -b;
            const bool insideLow = getExactSign(apexNumerator) == sa;
            const bool insideHigh = sa > 0 ? apexNumerator < a : apexNumerator > a;
            if (!insideLow || !insideHigh)
            {
                return false;
            }
            if (sd == 0)
            {
                isTangent = true;
                return false;
            }
            // f(t*) = -discriminant / a: opposite of the side of the ends when the roots are real and distinct
            const int sideOfEnds = anEnd1 == ON_POSITIVE_SIDE ? 1 : -1;
            return (-sd * sa) == -sideOfEnds;
        }

    private:
        /** @brief The coefficients of h(epsilon) . v(epsilon), up to cOrder*/
        static std::vector<S> values(const P& aPlane, const P* aPerturbation, const MathSymbolicPoint<P>& aPoint)
        {
            std::vector<S> result(cOrder + 1, S(0));
            for (size_t k = 0; k <= cOrder; k++)
            {
                if (k < aPoint.size())
                {
                    result[k] = aPlane.dot(aPoint.c[k]);
                }
                if (aPerturbation != nullptr && k >= 1 && k - 1 < aPoint.size())
                {
                    result[k] = result[k] + aPerturbation->dot(aPoint.c[k - 1]);
                }
            }
            return result;
        }
    };
}
