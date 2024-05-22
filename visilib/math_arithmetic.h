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
#include <math.h>
#include <float.h>

#include "visilib_core.h"

namespace visilib
{
    template<class S>
    class MathVector3_;

    /** @brief Wrapper for elementary arithmetic operations, allowing to use exact arithmetic operations.*/

    template<class S>
    class MathArithmetic
    {
    public:

        /** @brief Tests if the two numbers x and y have the same sign.*/
        static bool isSameSign(double x, double y)
        {
            return std::signbit(x) == std::signbit(y);
        }

        /** @brief Compute the axis aligned bounding box of a set of points.*/
        static void getMinMax(const std::vector<MathVector3_<S>>& aPoints, MathVector3_<S>& aMin, MathVector3_<S>& aMax);
      
        /** @brief Compute the axis aligned bounding box of a set of points.*/
        static void getMinMax(const S* aPoints, size_t pointCount, MathVector3_<S>& aMin, MathVector3_<S>& aMax);

        /** @brief Compute the absolute value of the number a.*/
        static S getAbs(S a);

        /** @brief Compute the square root of a.*/
        static S getSqrt(S a);

        /** @brief Tolerance value.*/
        static S Tolerance();

        /** @brief Tolerance value when clipping a polygon.*/
        static S GuardBandClipping();

        /** @brief Test is a point is finite (not a NaN).*/
        static bool isFinite(S a);

        /** @brief Test is a vector is finite (not containing NaN).*/
        static bool isFinite(MathVector3_<S> a);
    };

    template<> inline double MathArithmetic<double>::Tolerance() { return  1e-12; }
    template<> inline float MathArithmetic<float>::Tolerance() { return 1e-6f; }
    template<> inline double MathArithmetic<double>::GuardBandClipping() { return 1e-12; };
    template<> inline float MathArithmetic<float>::GuardBandClipping() { return 1e-6f; };

    template<> inline double MathArithmetic<double>::getAbs(double s);
    template<> inline float MathArithmetic<float>::getAbs(float s);

#ifdef EXACT_ARITHMETIC
    template<> exact MathArithmetic<exact>::getAbs(exact s);
#endif

    template<> inline double MathArithmetic<double>::getSqrt(double s);
    template<> inline float MathArithmetic<float>::getSqrt(float s);

#ifdef EXACT_ARITHMETIC
    template<> exact MathArithmetic<exact>::getSqrt(exact s);
#endif

    template<> inline bool MathArithmetic<double>::isFinite(double s);

    template<> inline bool MathArithmetic<float>::isFinite(float s);
#ifdef EXACT_ARITHMETIC
    template<> inline bool MathArithmetic<exact>::isFinite(exact s);
#endif

#ifdef EXACT_ARITHMETIC
    template<> inline exact MathArithmetic<exact>::Tolerance() { return leda_real(1e-20); }
#endif

    template<>
    inline double MathArithmetic<double>::getAbs(double s)
    {
        return std::fabs(s);
    }

    template<>
    inline float MathArithmetic<float>::getAbs(float s)
    {
        return std::fabs(s);
    }

#ifdef EXACT_ARITHMETIC
    template<>
    inline exact MathArithmetic<exact>::getAbs(exact s)
    {
        return CGAL::abs(s);
    }
#endif

    template<>
    inline double MathArithmetic<double>::getSqrt(double s)
    {
        return sqrt(s);
    }

    template<>
    inline float MathArithmetic<float>::getSqrt(float s)
    {
        return sqrt(s);
    }

#ifdef EXACT_ARITHMETIC
    template<>
    inline exact MathArithmetic<exact>::getSqrt(exact s)
    {
        //return approximate_sqrt(s);0
        return CGAL::sqrt(s);
    }
#endif

    template<>
    inline bool MathArithmetic<double>::isFinite(double s)
    {
        return isfinite(s);
    }

    template<>
    inline bool MathArithmetic<float>::isFinite(float s)
    {
        return isfinite(s);
    }

#ifdef EXACT_ARITHMETIC
    template<>
    inline bool MathArithmetic<exact>::isFinite(exact s)
    {
        return CGAL::is_finite(s);
    }
#endif

    template<class S>
    inline bool MathArithmetic<S>::isFinite(MathVector3_<S> a)
    {
        return isFinite(a.x) && isFinite(a.y) && isFinite(a.z);
    }

    template<class S>
    inline void  MathArithmetic<S>::getMinMax(const S* aPoints, size_t aPointCount, MathVector3_<S>& aMin, MathVector3_<S>& aMax)
    {
        const MathVector3_<S>* p = reinterpret_cast<const MathVector3_<S>*>(aPoints);
        aMin = MathVector3_<S>( FLT_MAX,  FLT_MAX,  FLT_MAX);
        aMax = MathVector3_<S>(-FLT_MAX, -FLT_MAX, -FLT_MAX);
        
        for (size_t i = 0; i < aPointCount; i++)
        {
            if (p[i].x < aMin.x) aMin.x = p[i].x;
            if (p[i].y < aMin.y) aMin.y = p[i].y;
            if (p[i].z < aMin.z) aMin.z = p[i].z;

            if (p[i].x > aMax.x) aMax.x = p[i].x;
            if (p[i].y > aMax.y) aMax.y = p[i].y;
            if (p[i].z > aMax.z) aMax.z = p[i].z;
        }
    }
    
    template<class S>
    inline void  MathArithmetic<S>::getMinMax(const std::vector<MathVector3_<S>>& aPoints, MathVector3_<S>& aMin, MathVector3_<S>& aMax)
    {
        const S* p = static_cast<const S*>(&(aPoints[0].x));
        MathArithmetic<S>::getMinMax(p, aPoints.size(), aMin, aMax);
    }
}
