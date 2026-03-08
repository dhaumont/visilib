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
#include <limits>
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
        static void getMinMax(const std::vector<MathVector3_<S> >& aPoints, MathVector3_<S>& aMin, MathVector3_<S>& aMax);

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

        static inline int bitsCount();

        static inline int validDecimalsCount();
        /** @brief Test is a point is finite (not a NaN).*/
        static bool isFinite(S a);

        /** @brief Test is a vector is finite (not containing NaN).*/
        static bool isFinite(MathVector3_<S> a);

        static double to_double(S a);
        static S getRandom();
    };
    template <> inline int MathArithmetic<double>::validDecimalsCount()  { return std::numeric_limits<double>::digits10 + 1; }
    template <> inline int MathArithmetic<float>::validDecimalsCount()  { return std::numeric_limits<float>::digits10 + 1; }
#ifdef ENABLE_LEDA
    template <> inline int MathArithmetic<MathLedaReal>::validDecimalsCount()  { return 200; }
#endif
#ifdef ENABLE_MPFR
    template <> inline int MathArithmetic<MathMpfr>::validDecimalsCount()  { return 200; }
#endif
#ifdef ENABLE_GMP
    template <> inline int MathArithmetic<MathGmpFloat>::validDecimalsCount()  { return 200; }
#endif

    template<> inline double MathArithmetic<double>::Tolerance() { return  1e-11; }
    template<> inline float MathArithmetic<float>::Tolerance() { return 1e-6f; }
#ifdef ENABLE_LEDA
    template<> inline MathLedaReal MathArithmetic<MathLedaReal>::Tolerance() 
    {
       int basis(10.0);
       int precision(MathArithmetic<MathLedaReal>::validDecimalsCount()-1);
       return pow(basis, -precision);
    }
#endif
#ifdef ENABLE_GMP
    template<> inline MathGmpFloat MathArithmetic<MathGmpFloat>::Tolerance()
    {
       int basis(10.0);
       int precision(MathArithmetic<MathGmpFloat>::validDecimalsCount()-1);
       return pow(basis, -precision);
    }
    template<> inline MathGmpRational MathArithmetic<MathGmpRational>::Tolerance()
    {
      return  MathGmpRational(1e-200);
    }
#endif
#ifdef ENABLE_MPFR
    template<> inline MathMpfr MathArithmetic<MathMpfr>::Tolerance()
    {
       int basis(10.0);
       int precision(MathArithmetic<MathMpfr>::validDecimalsCount()-1);
       return pow(basis, -precision);
    }
#endif
    template<> inline double MathArithmetic<double>::GuardBandClipping() { return 1e-12; };
    template<> inline float MathArithmetic<float>::GuardBandClipping() { return 1e-6f; };

    template<> inline double MathArithmetic<double>::getAbs(double s);
    template<> inline float MathArithmetic<float>::getAbs(float s);

    template <> inline int MathArithmetic<double>::bitsCount()  { return 64; }
    template <> inline int MathArithmetic<float>::bitsCount()  { return 32; }
#ifdef ENABLE_LEDA
    template <> inline int MathArithmetic<MathLedaReal>::bitsCount()
    {
      return 0;
    }
#endif
#ifdef ENABLE_MPFR
    template <> inline int MathArithmetic<MathMpfr>::bitsCount()
    {
      return MathMpfr::digits2bits(MathArithmetic<MathMpfr>::validDecimalsCount());
    }
#endif
#ifdef ENABLE_GMP
    template <> inline int MathArithmetic<MathGmpFloat>::bitsCount()
    {
      return MathGmpFloat::digits2bits(MathArithmetic<MathGmpFloat>::validDecimalsCount());
    }
    template <> inline int MathArithmetic<MathGmpRational>::bitsCount()
    {
      return 0;
    }
#endif
#ifdef ENABLE_LEDA
    template<> MathLedaReal MathArithmetic<MathLedaReal>::getAbs(MathLedaReal s);
#endif
#ifdef ENABLE_GMP
    template<> inline MathGmpFloat MathArithmetic<MathGmpFloat>::getAbs(MathGmpFloat s);
    template<> inline MathGmpRational MathArithmetic<MathGmpRational>::getAbs(MathGmpRational s);
#endif
#ifdef ENABLE_MPFR
    template<> inline MathMpfr MathArithmetic<MathMpfr>::getAbs(MathMpfr  s);
#endif

    template<> inline double MathArithmetic<double>::getSqrt(double s);
    template<> inline float MathArithmetic<float>::getSqrt(float s);

#ifdef ENABLE_LEDA
    template<> MathLedaReal MathArithmetic<MathLedaReal>::getSqrt(MathLedaReal s);
#endif
#ifdef ENABLE_GMP
    template<> inline MathGmpFloat MathArithmetic<MathGmpFloat>::getSqrt(MathGmpFloat  s);
    template<> inline MathGmpRational MathArithmetic<MathGmpRational>::getSqrt(MathGmpRational  s);
#endif
#ifdef ENABLE_MPFR
    template<> inline MathMpfr MathArithmetic<MathMpfr>::getSqrt(MathMpfr s);
#endif

    template<> inline bool MathArithmetic<double>::isFinite(double s);

    template<> inline bool MathArithmetic<float>::isFinite(float s);
#ifdef ENABLE_LEDA
    template<> inline bool MathArithmetic<MathLedaReal>::isFinite(MathLedaReal s);
#endif
#ifdef ENABLE_GMP
    template<> inline bool MathArithmetic<MathGmpFloat>::isFinite(MathGmpFloat s);
    template<> inline bool MathArithmetic<MathGmpRational>::isFinite(MathGmpRational  s);
#endif
#ifdef ENABLE_MPFR
    template<> inline bool MathArithmetic<MathMpfr>::isFinite(MathMpfr s);
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

#ifdef ENABLE_LEDA
    template<>
    inline MathLedaReal MathArithmetic<MathLedaReal>::getAbs(MathLedaReal s)
    {
        return abs(s);
    }
#endif
#ifdef ENABLE_GMP
    template<>
    inline MathGmpFloat MathArithmetic<MathGmpFloat>::getAbs(MathGmpFloat  s)
    {
        return s.abs();
    }

    template<>
    inline MathGmpRational MathArithmetic<MathGmpRational>::getAbs(MathGmpRational s)
    {
        return s.abs();
    }
#endif
#ifdef ENABLE_MPFR
    template<>
    inline MathMpfr MathArithmetic<MathMpfr>::getAbs(MathMpfr s)
    {
        return s.abs();
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

#ifdef ENABLE_LEDA
    template<>
    inline MathLedaReal MathArithmetic<MathLedaReal>::getSqrt(MathLedaReal s)
    {
        //return approximate_sqrt(s);0
        return sqrt(s);
    }
#endif

#ifdef ENABLE_GMP
    template<>
    inline MathGmpFloat MathArithmetic<MathGmpFloat>::getSqrt(MathGmpFloat s)
    {
        return s.sqrt();
    }
    template<>
    inline MathGmpRational MathArithmetic<MathGmpRational>::getSqrt(MathGmpRational  s)
    {
        return s.sqrt();
    }
#endif

#ifdef ENABLE_MPFR
    template<>
    inline MathMpfr MathArithmetic<MathMpfr>::getSqrt(MathMpfr s)
    {
        return s.sqrt();
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

#ifdef ENABLE_LEDA
    template<>
    inline bool MathArithmetic<MathLedaReal>::isFinite(MathLedaReal s)
    {
        return true;
    }
#endif
#ifdef ENABLE_GMP
    template<>
    inline bool MathArithmetic<MathGmpFloat>::isFinite(MathGmpFloat s)
    {
        return s.isfinite();
    }
    template<>
    inline bool MathArithmetic<MathGmpRational>::isFinite(MathGmpRational s)
    {
        return s.isfinite();
    }
#endif
#ifdef ENABLE_MPFR
    template<>
    inline bool MathArithmetic<MathMpfr>::isFinite(MathMpfr s)
    {
        return s.isfinite();
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

    template<>
    inline double MathArithmetic<float>::to_double(float v)
    {
        return double(v);
    }

    template<>
    inline double MathArithmetic<double>::to_double(double v)
    {
        return v;
    }
#ifdef ENABLE_LEDA
    template<>
    inline double MathArithmetic<MathLedaReal>::to_double(MathLedaReal v)
    {
        return v.to_double();
    }
#endif
#ifdef ENABLE_MPFR
    template<>
    inline double MathArithmetic<MathMpfr>::to_double(MathMpfr v)
    {
        return v.to_double();
    }
#endif
#ifdef ENABLE_GMP
    template<>
    inline double MathArithmetic<MathGmpFloat>::to_double(MathGmpFloat v)
    {
        return v.to_double();
    }
    template<>
    inline double MathArithmetic<MathGmpRational>::to_double(MathGmpRational v)
    {
        return v.to_double();
    }
#endif
    template<class S>
    inline S MathArithmetic<S>::getRandom()
    {
        return (S)rand() / (S)(RAND_MAX);
    }

}
