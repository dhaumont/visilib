/*
Visilib, an open source library for exact visibility computation.
Copyright(C) 2020 by Denis Haumont

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

#include <cmath>
#include <functional>
#include <iostream>

namespace visilib
{
    template<class S>
    class MathArithmetic;

    template<class S>
    class MathVector2_;

    template<class S>
    class MathPlane3_;

    template<class S>
    class MathVector3_;

    typedef MathVector2_<float> MathVector2f;
    typedef MathVector2_<double> MathVector2d;
    typedef MathVector2_<size_t> MathVector2i;


    /** @brief Represents a vector in 2D euclidian space.
       @param S the coordinates representation (float, double, exact)
    */

    template<class S>
    class MathVector2_
    {
    public:
        MathVector2_() : x(0), y(0) {}
        MathVector2_(const MathVector2_& aV) : x(aV.x), y(aV.y) {}
        MathVector2_(const MathVector2_& aV1, const MathVector2_& aV2) : x(aV1.x - aV2.x), y(aV1.y - aV2.y) {}
        MathVector2_(S aX, S aY) : x(aX), y(aY) {}

        MathVector2_& operator =  (const MathVector2_& v) { x = v.x;  y = v.y;  return *this; }
        MathVector2_& operator += (const MathVector2_& v) { x += v.x; y += v.y; return *this; }
        MathVector2_& operator -= (const MathVector2_& v) { x -= v.x; y -= v.y; return *this; }
        MathVector2_& operator *= (const MathVector2_& v) { x *= v.x; y *= v.y; return *this; }
        MathVector2_& operator /= (const MathVector2_& v) { x /= v.x; y /= v.y; return *this; }
        MathVector2_& operator *= (S f) { x *= f;  y *= f; return *this; }
        MathVector2_& operator /= (S f) { x /= f;  y /= f; return *this; }

        MathVector2_ operator * (S f) { return MathVector2_(x * f, y * f); }
        MathVector2_ operator -(const MathVector2_& v) const { return MathVector2_(x - v.x, y - v.y); }
        MathVector2_ operator +(const MathVector2_ & v) const { return MathVector2_(x + v.x, y + v.y); }

        bool operator==(const MathVector2_ & v) const { return(x == v.x && y == v.y); }
        bool operator!=(const MathVector2_ & v) const { return(x != v.x || y != v.y); }

        S* operator &() { return (&x); }
        const S* operator &() const { return (&x); }

        S dot(const MathVector2_ & a)const { return a.x* x + a.y * y; }

        bool isZero(S anEps)const
        {
            return (std::fabs(x) <= anEps && std::fabs(y) <= anEps);
        }

        S getSquaredNorm()const { return (dot(*this)); }
        S getNorm()const { return MathArithmetic<S>::getSqrt(getSquaredNorm()); }

        friend std::ostream& operator<<(std::ostream & os, const MathVector2_ & v)
        {
            os << "(" << v.x << "," << v.y << ")";
            return os;
        }

        static MathVector3_<S> cross(const MathVector2_ & v1, const MathVector2_ & v2);

        S normalize()
        {
            S d = getNorm();

#ifdef _DEBUG
            // V_ASSERT (std::fabs(d)>MathArithmetic::Tolerance());
#endif

            if (std::fabs(d) > MathArithmetic<S>::Tolerance())
            {
                S inv = (S)1.0 / d;
                x *= inv;
                y *= inv;
            }

            return d;
        }
        MathVector2_ getNormalized() const
        {
            MathVector2_ result(*this);
            result.normalize();
            return result;
        }

        template<typename F>
        MathVector2_<F> to() const;

        static MathVector2_ Zero() { return MathVector2_((S)0.0, (S)0.0); }
        static MathVector2_ X() { return MathVector2_((S)1.0, (S)0.0); }
        static MathVector2_ Y() { return MathVector2_((S)0.0, (S)1.0); }

        void interp(const MathVector2_<S> & myV1, const MathVector2_<S> & myV2, const MathPlane3_<S> & aPlane);

        /**@brief A hash associated to an edge stored in a MathVector2i*/
        struct Hash
        {
        public:
            /**@brief Compute the haskey of the edge*/
            size_t operator()(const MathVector2i& a) const
            {
                return (size_t)std::hash<S>()(a.x) ^ (size_t)std::hash<S>()(a.y);
            }
        };

        S x, y;
    };

    template<class S>
    MathVector3_<S> MathVector2_<S>::cross(const MathVector2_ & v1, const MathVector2_ & v2)
    {
        return MathVector3_<S>(0, 0, v1.x * v2.y - v1.y * v2.x);
    }
}
