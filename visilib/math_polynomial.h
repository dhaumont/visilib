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
#include <vector>

#include "math_arithmetic.h"

namespace visilib
{
    template<class S>
    class MathPolynomial_;

    typedef MathPolynomial_<double> MathPolynomiald;

    /** @brief A polynomial of one variable, c_0 + c_1 t + ... + c_n t^n, and the real roots it has on an interval.
        @param S the coordinates representation (float, double, exact)*/
    template<class S>
    class MathPolynomial_
    {
    public:
        MathPolynomial_()
        {
        }

        /** @brief The polynomial of the given coefficients, the constant term first.*/
        explicit MathPolynomial_(const std::vector<S>& someCoefficients)
            : mCoefficients(someCoefficients)
        {
        }

        S evaluate(const S& t) const
        {
            S result = S(0);
            for (size_t i = mCoefficients.size(); i > 0; i--)
            {
                result = result * t + mCoefficients[i - 1];
            }
            return result;
        }

        MathPolynomial_ getDerivative() const
        {
            std::vector<S> d;
            for (size_t i = 1; i < mCoefficients.size(); i++)
            {
                d.push_back(mCoefficients[i] * S(i));
            }
            return MathPolynomial_(d);
        }

        /** @brief The greatest of the absolute values of the coefficients.*/
        S getScale() const
        {
            S scale = S(0);
            for (const S& c : mCoefficients)
            {
                scale = std::max(scale, MathArithmetic<S>::getAbs(c));
            }
            return scale;
        }

        /** @brief The degree once the leading coefficients below aRelativeTolerance times the scale are dropped (0 for a null polynomial).*/
        size_t getDegree(const S& aRelativeTolerance) const
        {
            const S threshold = getScale() * aRelativeTolerance;
            size_t degree = mCoefficients.size();
            while (degree > 0 && MathArithmetic<S>::getAbs(mCoefficients[degree - 1]) <= threshold)
            {
                degree--;
            }
            return degree == 0 ? 0 : degree - 1;
        }

        /** @brief Whether every coefficient is below aRelativeTolerance times aReference.*/
        bool isNull(const S& aReference, const S& aRelativeTolerance) const
        {
            return getScale() <= aReference * aRelativeTolerance;
        }

        const std::vector<S>& getCoefficients() const { return mCoefficients; }

        /** @brief Append the real roots of the polynomial that lie in [aLow, aHigh] to roots, sorted, each once.

        The coefficients below aRelativeTolerance times the scale are dropped from the top. Up to the second degree the
        roots are those of the closed form; beyond, the roots of the derivative, found the same way, split the interval in
        pieces where the polynomial is monotone, and each piece with a change of sign holds one root, found by bisection.
        A root of even multiplicity, where the sign does not change, is found only as a root of the derivative where the
        value vanishes. A null polynomial has no isolated root and appends nothing.*/
        void findRoots(const S& aLow, const S& aHigh, const S& aRelativeTolerance, std::vector<S>& roots) const
        {
            std::vector<S> found;
            findRootsInternal(aLow, aHigh, aRelativeTolerance, found);
            std::sort(found.begin(), found.end());
            const S merge = (aHigh - aLow) * aRelativeTolerance;
            for (size_t i = 0; i < found.size(); i++)
            {
                if (i == 0 || found[i] - found[i - 1] > merge)
                {
                    roots.push_back(found[i]);
                }
            }
        }

    private:
        void findRootsInternal(const S& aLow, const S& aHigh, const S& aRelativeTolerance, std::vector<S>& roots) const
        {
            const size_t degree = getDegree(aRelativeTolerance);
            if (degree == 0)
            {
                return;
            }

            const S scale = getScale();
            if (degree <= 2)
            {
                // a t^2 + 2 b t + c with b the half of the linear coefficient
                const S a = degree == 2 ? mCoefficients[2] / scale : S(0);
                const S b = mCoefficients[1] / scale / S(2);
                const S c = mCoefficients[0] / scale;
                std::vector<S> candidates;
                MathArithmetic<S>::solveHalfQuadratic(a, b, c, aRelativeTolerance, candidates);
                for (const S& t : candidates)
                {
                    if (t >= aLow && t <= aHigh)
                    {
                        roots.push_back(t);
                    }
                }
                return;
            }

            // the pieces where the polynomial is monotone, bounded by the roots of the derivative
            std::vector<S> bounds;
            bounds.push_back(aLow);
            getDerivative().findRoots(aLow, aHigh, aRelativeTolerance, bounds);
            bounds.push_back(aHigh);

            const S zero = scale * aRelativeTolerance;
            for (size_t i = 0; i + 1 < bounds.size(); i++)
            {
                S lo = bounds[i];
                S hi = bounds[i + 1];
                S flo = evaluate(lo);
                S fhi = evaluate(hi);

                if (MathArithmetic<S>::getAbs(flo) <= zero)
                {
                    roots.push_back(lo); // a root at a critical point, possibly of even multiplicity
                    continue;
                }
                if (i + 2 == bounds.size() && MathArithmetic<S>::getAbs(fhi) <= zero)
                {
                    roots.push_back(hi);
                    continue;
                }
                if ((flo < S(0)) == (fhi < S(0)))
                {
                    continue;
                }
                for (int iteration = 0; iteration < 200 && hi - lo > (aHigh - aLow) * std::numeric_limits<double>::epsilon(); iteration++)
                {
                    const S mid = (lo + hi) / S(2);
                    const S fmid = evaluate(mid);
                    if ((fmid < S(0)) == (flo < S(0)))
                    {
                        lo = mid;
                        flo = fmid;
                    }
                    else
                    {
                        hi = mid;
                    }
                }
                roots.push_back((lo + hi) / S(2));
            }
        }

        std::vector<S> mCoefficients; /**< @brief The coefficients, the constant term first*/
    };
}
