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

#include <cmath>
#include <cstdint>

#include "math_vector_3.h"

namespace visilib
{
    /** @brief A signed integer of 128 bits, enough for the products of three coordinates of the grid (see MathIntegerGrid).

    Only what the predicates need: the product of two 64 bit integers, the sum, the difference and the sign. The
    compilers that have a native 128 bit integer use it; the others get a portable version on two 64 bit words.*/
    class MathInt128
    {
    public:
        MathInt128() : MathInt128(int64_t(0)) {}

#ifdef __SIZEOF_INT128__
        explicit MathInt128(int64_t v) : mValue(v) {}

        static MathInt128 multiply(int64_t a, int64_t b)
        {
            MathInt128 r;
            r.mValue = (__int128)a * (__int128)b;
            return r;
        }
        /** @brief The product of a 128 bit integer by a 64 bit one: the caller guarantees it fits 128 bits.*/
        static MathInt128 multiply(const MathInt128& a, int64_t b)
        {
            MathInt128 r;
            r.mValue = a.mValue * (__int128)b;
            return r;
        }
        MathInt128 operator+(const MathInt128& o) const { MathInt128 r; r.mValue = mValue + o.mValue; return r; }
        MathInt128 operator-(const MathInt128& o) const { MathInt128 r; r.mValue = mValue - o.mValue; return r; }
        int sign() const { return mValue > 0 ? 1 : (mValue < 0 ? -1 : 0); }

    private:
        __int128 mValue;
#else
        explicit MathInt128(int64_t v) : mLow(uint64_t(v)), mHigh(v < 0 ? ~uint64_t(0) : 0) {}

        static MathInt128 multiply(int64_t a, int64_t b)
        {
            const bool negative = (a < 0) != (b < 0);
            MathInt128 r = unsignedProduct(magnitude(a), magnitude(b));
            return negative ? r.negate() : r;
        }

        /** @brief The product of a 128 bit integer by a 64 bit one: the caller guarantees it fits 128 bits.*/
        static MathInt128 multiply(const MathInt128& a, int64_t b)
        {
            const bool negativeA = a.sign() < 0;
            const MathInt128 ua = negativeA ? a.negate() : a;
            const uint64_t ub = magnitude(b);
            // (high 2^64 + low) b, of which only the low 128 bits are kept, the product fitting them
            MathInt128 r = unsignedProduct(ua.mLow, ub);
            r.mHigh += ua.mHigh * ub;
            return (negativeA != (b < 0)) ? r.negate() : r;
        }
        MathInt128 operator+(const MathInt128& o) const
        {
            MathInt128 r;
            r.mLow = mLow + o.mLow;
            r.mHigh = mHigh + o.mHigh + (r.mLow < mLow ? 1 : 0);
            return r;
        }
        MathInt128 operator-(const MathInt128& o) const { return *this + o.negate(); }
        int sign() const
        {
            if (int64_t(mHigh) < 0) return -1;
            return (mHigh != 0 || mLow != 0) ? 1 : 0;
        }

    private:
        static uint64_t magnitude(int64_t v)
        {
            return v < 0 ? uint64_t(0) - uint64_t(v) : uint64_t(v);
        }

        // the product of two unsigned 64 bit integers on 32 bit halves
        static MathInt128 unsignedProduct(uint64_t ua, uint64_t ub)
        {
            const uint64_t a0 = ua & 0xffffffffu, a1 = ua >> 32;
            const uint64_t b0 = ub & 0xffffffffu, b1 = ub >> 32;
            const uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
            const uint64_t middle = (p00 >> 32) + (p01 & 0xffffffffu) + (p10 & 0xffffffffu);
            MathInt128 r;
            r.mLow = (p00 & 0xffffffffu) | (middle << 32);
            r.mHigh = p11 + (p01 >> 32) + (p10 >> 32) + (middle >> 32);
            return r;
        }

        MathInt128 negate() const
        {
            MathInt128 r;
            r.mLow = ~mLow + 1;
            r.mHigh = ~mHigh + (r.mLow == 0 ? 1 : 0);
            return r;
        }
        uint64_t mLow;
        uint64_t mHigh;
#endif
    };

    /** @brief A point of the integer grid.*/
    struct MathIntegerPoint
    {
        int64_t x, y, z;

        bool operator==(const MathIntegerPoint& o) const { return x == o.x && y == o.y && z == o.z; }
        bool operator!=(const MathIntegerPoint& o) const { return !(*this == o); }
    };

    /** @brief A regular grid of the space, the coordinates of the input are snapped to, so that the predicates on them are
    exact with fixed width integers (the plane based integer arithmetic of Nehring-Wirxel, Trettner and Kobbelt).

    The step is a power of two and the coordinates have at most cBits bits. The grid is fine enough to leave the floats
    of the input as they are: a float of magnitude at least 2^24 steps is already a multiple of the step, and only the
    smaller ones are rounded, by half a step at most, to values that are floats too. The collinearities and coplanarities
    of the input survive the snapping, which a coarse grid would turn into tiny kinks. The difference of two coordinates
    takes cBits + 1 bits; the orientation of four points is a determinant of such differences, products of three factors,
    3 (cBits + 1) bits, plus three bits for the sums: 126 bits, which a 128 bit integer holds.

    The grid is centered on the origin and holds the cube [-extent, extent]^3: a coordinate outside cannot be snapped,
    and the points built from it have no exact predicate (see isOnGrid).*/
    class MathIntegerGrid
    {
    public:
        static constexpr int cBits = 40;

        MathIntegerGrid() : mStep(0.0), mValid(false) {}

        /** @brief The grid holding the cube [-anExtent, anExtent]^3, with the finest power of two step that fits cBits bits.*/
        explicit MathIntegerGrid(double anExtent)
        {
            mValid = anExtent > 0.0 && std::isfinite(anExtent);
            mStep = mValid ? std::ldexp(1.0, int(std::ceil(std::log2(anExtent))) - (cBits - 1)) : 0.0;
        }

        bool isValid() const { return mValid; }
        double getStep() const { return mStep; }
        double getExtent() const { return std::ldexp(mStep, cBits - 1); }

        /** @brief The coordinate snapped to the nearest grid value (unchanged when the grid is not valid).*/
        float snap(float aValue) const
        {
            if (!mValid)
            {
                return aValue;
            }
            return float(std::nearbyint(double(aValue) / mStep) * mStep);
        }

        MathVector3f snap(const MathVector3f& p) const
        {
            return MathVector3f(snap(p.x), snap(p.y), snap(p.z));
        }

        /** @brief Whether the point lies exactly on the grid, inside its extent: only such points have exact predicates.*/
        template<class V>
        bool isOnGrid(const V& p) const
        {
            if (!mValid)
            {
                return false;
            }
            const double c[3] = { double(p.x), double(p.y), double(p.z) };
            const double limit = std::ldexp(1.0, cBits - 1);
            for (int i = 0; i < 3; i++)
            {
                const double q = c[i] / mStep;
                if (q != std::nearbyint(q) || std::fabs(q) > limit)
                {
                    return false;
                }
            }
            return true;
        }

        /** @brief The integer coordinates of a point of the grid (see isOnGrid).*/
        template<class V>
        MathIntegerPoint toInteger(const V& p) const
        {
            return MathIntegerPoint{ int64_t(std::nearbyint(double(p.x) / mStep)),
                                     int64_t(std::nearbyint(double(p.y) / mStep)),
                                     int64_t(std::nearbyint(double(p.z) / mStep)) };
        }

    private:
        double mStep;
        bool mValid;
    };

    /** @brief Exact predicates on the points of an integer grid.*/
    class MathIntegerPredicates
    {
    public:
        /** @brief The sign of the orientation of four points: the sign of det(b - a, c - a, d - a). Zero when coplanar.*/
        static int orient3d(const MathIntegerPoint& a, const MathIntegerPoint& b, const MathIntegerPoint& c, const MathIntegerPoint& d)
        {
            const int64_t bx = b.x - a.x, by = b.y - a.y, bz = b.z - a.z;
            const int64_t cx = c.x - a.x, cy = c.y - a.y, cz = c.z - a.z;
            const int64_t dx = d.x - a.x, dy = d.y - a.y, dz = d.z - a.z;

            // cofactors of the first row: products of two differences, 2 (cBits + 1) + 1 bits, in 128 bits
            const MathInt128 m0 = MathInt128::multiply(cy, dz) - MathInt128::multiply(cz, dy);
            const MathInt128 m1 = MathInt128::multiply(cz, dx) - MathInt128::multiply(cx, dz);
            const MathInt128 m2 = MathInt128::multiply(cx, dy) - MathInt128::multiply(cy, dx);

            const MathInt128 det = MathInt128::multiply(m0, bx) + MathInt128::multiply(m1, by) + MathInt128::multiply(m2, bz);
            return det.sign();
        }

        /** @brief The sign of the Klein product of the oriented lines a->b and c->d: the sign of the Plucker side product
        of (b - a, a x b) and (d - c, c x d), as MathPlucker6::dot computes it. It is the orientation of the four points,
        and vanishes exactly when the two lines are coplanar (meeting or parallel).*/
        static int kleinSign(const MathIntegerPoint& a, const MathIntegerPoint& b, const MathIntegerPoint& c, const MathIntegerPoint& d)
        {
            return orient3d(a, b, c, d);
        }

        /** @brief Whether the three points are on one line: the cross product of b - a and c - a vanishes.*/
        static bool areCollinear(const MathIntegerPoint& a, const MathIntegerPoint& b, const MathIntegerPoint& c)
        {
            const int64_t bx = b.x - a.x, by = b.y - a.y, bz = b.z - a.z;
            const int64_t cx = c.x - a.x, cy = c.y - a.y, cz = c.z - a.z;
            return (MathInt128::multiply(by, cz) - MathInt128::multiply(bz, cy)).sign() == 0
                && (MathInt128::multiply(bz, cx) - MathInt128::multiply(bx, cz)).sign() == 0
                && (MathInt128::multiply(bx, cy) - MathInt128::multiply(by, cx)).sign() == 0;
        }
    };
}
