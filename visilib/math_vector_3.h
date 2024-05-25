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

#include "math_arithmetic.h"
#include <iostream>

namespace visilib
{
    template<class S>
    class MathVector3_;

    template<class S>
    class MathPlane3_;

    //typedef MathVector3_<S> MathVector3_<S>;
    typedef MathVector3_<float> MathVector3f;
    typedef MathVector3_<double> MathVector3d;


    /** @brief Represents a vector in 3D euclidian space.
       @param S the coordinates representation (float, double, exact)
    */

    template<class S>
    class MathVector3_
    {
    public:
        MathVector3_() : x(0), y(0), z(0) {}
        MathVector3_(const MathVector3_& aV) : x(aV.x), y(aV.y), z(aV.z) {}
        MathVector3_(const MathVector3_& aV1, const MathVector3_& aV2) : x(aV1.x - aV2.x), y(aV1.y - aV2.y), z(aV1.z - aV2.z) {}
        MathVector3_(S aX, S aY, S aZ) : x(aX), y(aY), z(aZ) {}
        MathVector3_(const S* a) : x(a[0]), y(a[1]), z(a[2]) {}

        MathVector3_& operator =  (const MathVector3_& v) { x = v.x;  y = v.y;  z = v.z;  return *this; }
        MathVector3_& operator += (const MathVector3_& v) { x += v.x; y += v.y; z += v.z; return *this; }
        MathVector3_& operator -= (const MathVector3_& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
        MathVector3_& operator *= (const MathVector3_& v) { x *= v.x; y *= v.y; z *= v.z; return *this; }
        MathVector3_& operator /= (const MathVector3_& v) { x /= v.x; y /= v.y; z /= v.z; return *this; }
        MathVector3_& operator *= (S f) { x *= f;  y *= f; z *= f; return *this; }
        MathVector3_& operator /= (S f) { x /= f;  y /= f; z /= f; return *this; }

        MathVector3_ operator * (S f) { return MathVector3_(x * f, y * f, z * f); }
        MathVector3_ operator -(const MathVector3_& v) const { return MathVector3_(x - v.x, y - v.y, z - v.z); }
        MathVector3_ operator +(const MathVector3_ & v) const { return MathVector3_(x + v.x, y + v.y, z + v.z); }

        bool operator==(const MathVector3_ & v) const { return(x == v.x && y == v.y && z == v.z); }
        bool operator!=(const MathVector3_ & v) const { return(x != v.x || y != v.y || z != v.z); }

        S* operator &() { return (&x); }
        const S* operator &() const { return (&x); }

        S operator[] (int i) const
        {
            return reinterpret_cast<const S*>(&x)[i];
        }

        friend std::ostream& operator<<(std::ostream & os, const MathVector3_ & v)
        {
            os << "(" << v.x << "," << v.y << "," << v.z << ")";
            return os;
        }

        S dot(const MathVector3_ & a)const { return a.x* x + a.y * y + a.z * z; }

        bool isZero(S anEps)const;

        S getSquaredNorm()const { return (dot(*this)); }
        S getNorm()const;

        static MathVector3_ cross(const MathVector3_ & v1, const MathVector3_ & v2)
        {
            return MathVector3_(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
        }

        S normalize();

        MathVector3_ getNormalized() const
        {
            MathVector3_ result(*this);
            result.normalize();
            return result;
        }

        static MathVector3_ Zero() { return MathVector3_((S)0.0, (S)0.0, (S)0.0); }
        static MathVector3_ X() { return MathVector3_((S)1.0, (S)0.0, (S)0.0); }
        static MathVector3_ Y() { return MathVector3_((S)0.0, (S)1.0, (S)0.0); }
        static MathVector3_ Z() { return MathVector3_((S)0.0, (S)0.0, (S)1.0); }

        S x, y, z;
    };

    template<typename To, typename From> To convert(const From & f);

    template<>
    inline MathVector3f convert(const MathVector3d & f)
    {
        return MathVector3f((float)f.x, (float)f.y, (float)f.z);
    }

    template<>
    inline MathVector3d convert(const MathVector3f & f)
    {
        return MathVector3d((double)f.x, (double)f.y, (double)f.z);
    }

    template<>
    inline MathVector3f convert(const MathVector3f & f)
    {
        return MathVector3f(f.x, f.y, f.z);
    }

    template<>
    inline MathVector3d convert(const MathVector3d & f)
    {
        return MathVector3d(f.x, f.y, f.z);
    }

#ifdef EXACT_ARITHMETIC
    template<>
    inline MathVector3d convert(const MathVector3_<exact> & f)
    {
        return MathVector3d(to_double(f.x), to_double(f.y), to_double(f.z));
    }

    template<>
    inline MathVector3_<exact> convert(const MathVector3d & f)
    {
        return MathVector3_<exact>(f.x, f.y, f.z);
    }
#endif

    template<class S>
    inline bool MathVector3_<S>::isZero(S anEps)const
    {
        return (MathArithmetic<S>::getAbs(x) <= anEps && MathArithmetic<S>::getAbs(y) <= anEps && MathArithmetic<S>::getAbs(z) <= anEps);
    }

    template<class S>
    inline S MathVector3_<S>::normalize()
    {
        S d = getNorm();

#ifdef _DEBUG
        // V_ASSERT (std::fabs(d)>MathArithmetic::Tolerance());
#endif

        if (MathArithmetic<S>::getAbs(d) > MathArithmetic<S>::Tolerance())
        {
            S inv = (S)1.0 / d;
            x *= inv;
            y *= inv;
            z *= inv;
        }
        else 
        {
             x = y = z = 0;
        }
        return d;
    }

    template<class S>
    S MathVector3_<S>::getNorm()const
    {
        return MathArithmetic<S>::getSqrt(getSquaredNorm());
    }
}
