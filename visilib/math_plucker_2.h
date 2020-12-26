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

#include "math_vector_2.h"
#include "visilib_core.h"

namespace visilib
{

    /** @brief Represents a point in 2D Plucker space.

    The points in Plucker space are encoded as two 2D vectors.

       @param S the coordinates representation (float, double, exact)
    */
    template<class S>
    class MathPlucker2
    {
    public:

        MathPlucker2()
        {
        }

        MathPlucker2(S x0, S x1, S x2)
            :mDirection(x0, x1),
            mLocation(x2)
        {
        }

        MathPlucker2(const MathPlucker2& v)
        {
            mDirection = v.mDirection;
            mLocation = v.mLocation;
        }

        MathPlucker2(const MathVector2d& aBegin, const MathVector2d& anEnd)
        {
            MathVector2_<S> myBegin = aBegin.to<S>();
            MathVector2_<S> myEnd = anEnd.to<S>();

            mDirection = myEnd - myBegin;
            mLocation = MathVector2_<S>::cross(myBegin, myEnd).x;
        }

        MathPlucker2(const MathVector3d& aBegin, const MathVector3d& anEnd)
        {
            MathVector2_<S> myBegin = MathVector2_<S>(aBegin.x, aBegin.y);
            MathVector2_<S> myEnd = MathVector2_<S>(anEnd.x, anEnd.y);

            mDirection = myEnd - myBegin;
            mLocation = MathVector2_<S>::cross(myBegin, myEnd).x;
        }

        MathPlucker2& operator -= (const MathPlucker2& v)
        {
            mLocation -= v.mLocation;
            mDirection -= v.mDirection;

            return *this;
        }

        MathPlucker2& operator += (const MathPlucker2& v)
        {
            mLocation += v.mLocation;
            mDirection += v.mDirection;

            return *this;
        }

        MathPlucker2& operator *= (S aNumber)
        {
            mLocation *= aNumber;
            mDirection *= aNumber;

            return *this;
        }

        MathPlucker2 operator +(const MathPlucker2& v) const
        {
            MathPlucker2 result(*this);
            result += v;

            return result;
        }

        MathPlucker2 operator -(const MathPlucker2& v) const
        {
            MathPlucker2 result(*this);
            result -= v;

            return result;
        }

        MathPlucker2 operator / (S aNumber) const
        {
            MathPlucker2 result(*this);
            result.mLocation /= aNumber;
            result.mDirection /= aNumber;

            return result;
        }

        MathPlucker2 operator * (S aNumber) const
        {
            MathPlucker2 result(*this);
            result.mLocation *= aNumber;
            result.mDirection *= aNumber;

            return result;
        }

        const MathVector2_<S>& getDirection()const
        {
            return mDirection;
        }

        S getLocation()const
        {
            return mLocation;
        }

        MathPlucker2 getNormalized() const
        {
            MathPlucker2 result(*this);
            S myDiv = mDirection.getNorm();
            V_ASSERT(fabs(myDiv) > MathArithmetic<S>::Tolerance());

            result.mDirection /= myDiv;
            result.mLocation /= myDiv;
            return result;
        }

        S dot(const MathPlucker2<S>& p2) const
        {
            return mDirection.dot(p2.mDirection) + mLocation;
        }

        static MathPlucker2<S>  Zero()
        {
            return MathPlucker2(0.0, 0.0, 0.0);
        }
    private:
        MathVector2_<S> mDirection;
        S mLocation;
    };

    template<class S>
    std::ostream& operator<<(std::ostream& os, const MathPlucker2<S>& v)
    {
        os << v.mDirection << "," << v.mLocation;
        return os;
    }
}
