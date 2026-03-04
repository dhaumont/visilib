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

#include "visilib_core.h"
#include "math_arithmetic.h"
#include "math_vector_3.h"
#include "math_plane_3.h"

namespace visilib
{

    /** @brief Represents a point in 6D Plucker space.

    The 6D points in Plucker space are encoded as two 3D vectors.
        @param S: the coordinates representation (float, double, exact)
    */
    template<class S>
    class MathPlucker6
    {
    public:

        MathPlucker6()
        {
        }

        /**@brief Construct a plucker point from its 6 coefficients in Plucker space*/
        MathPlucker6(S x0, S x1, S x2, S x3, S x4, S x5)
            :mDirection(x0, x1, x2),
            mLocation(x3, x4, x5)
        {
            V_ASSERT(MathArithmetic<S>::isFinite(getDirection()));
            V_ASSERT(MathArithmetic<S>::isFinite(getLocation()));
        }

        /**@brief Construct a plucker point from another point*/
        MathPlucker6(const MathPlucker6<S>& v)
        {
            mDirection = v.mDirection;
            mLocation = v.mLocation;

            V_ASSERT(MathArithmetic<S>::isFinite(getDirection()));
            V_ASSERT(MathArithmetic<S>::isFinite(getLocation()));
        }

        /**@brief Construct a plucker point from a 3D line defined by two points*/
        MathPlucker6(const MathVector3d& aBegin, const MathVector3d& anEnd)
        {
            MathVector3_<S> myBegin = convert<MathVector3_<S> >(aBegin);
            MathVector3_<S> myEnd = convert<MathVector3_<S> >(anEnd);

            mDirection = myEnd - myBegin;
            mLocation = MathVector3_<S>::cross(myBegin, myEnd);

            V_ASSERT(MathArithmetic<S>::isFinite(getDirection()));
            V_ASSERT(MathArithmetic<S>::isFinite(getLocation()));
        }

        /**@brief Construct a plucker point from a 3D line defined by two points*/
        MathPlucker6(const MathVector3f& aBegin, const MathVector3f& anEnd)
        {
            MathVector3_<S> myBegin = convert<MathVector3_<S>>(aBegin);
            MathVector3_<S> myEnd = convert<MathVector3_<S>>(anEnd);

            mDirection = myEnd - myBegin;
            mLocation = MathVector3_<S>::cross(myBegin, myEnd);

            V_ASSERT(MathArithmetic<S>::isFinite(getDirection()));
            V_ASSERT(MathArithmetic<S>::isFinite(getLocation()));
        }

        /**@brief Subtract another point in Plucker space*/
        MathPlucker6<S>& operator -= (const MathPlucker6<S>& v)
        {
            mLocation -= v.mLocation;
            mDirection -= v.mDirection;

            return *this;
        }

        /**@brief Add another point in Plucker space*/
        MathPlucker6<S>& operator += (const MathPlucker6<S>& v)
        {
            mLocation += v.mLocation;
            mDirection += v.mDirection;

            return *this;
        }

        /**@brief Multiply all the coordinates by a provided number in Plucker space*/
        MathPlucker6<S>& operator *= (S aNumber)
        {
            mLocation *= aNumber;
            mDirection *= aNumber;

            return *this;
        }

        /**@brief Return the addition of two points in Plucker space*/
        MathPlucker6<S>  operator +(const MathPlucker6<S>& v) const
        {
            MathPlucker6<S>  result(*this);
            result += v;

            return result;
        }

        /**@brief Return the subtraction of two points in Plucker space*/
        MathPlucker6<S>  operator -(const MathPlucker6<S>& v) const
        {
            MathPlucker6<S>  result(*this);
            result -= v;

            return result;
        }

        /**@brief Return a point with its coordinates divided by the given number in Plucker space*/
        MathPlucker6<S>  operator / (S aNumber) const
        {
            MathPlucker6<S>  result(*this);
            result.mLocation /= aNumber;
            result.mDirection /= aNumber;

            return result;
        }

        /**@brief Rerturn a point with its coordinates multiplied by the given number in Plucker space*/
        MathPlucker6<S>  operator * (S aNumber) const
        {
            MathPlucker6<S>  result(*this);
            result.mLocation *= aNumber;
            result.mDirection *= aNumber;

            return result;
        }

        /**@brief Return the directional part of the Plucker point*/
        const MathVector3_<S>& getDirection()const
        {
            return mDirection;
        }


        /**@brief Return the location part of the Plucker point*/
        const MathVector3_<S>& getLocation()const
        {
            return mLocation;
        }

        /** @brief Compute the Plucker dot product of the point with another point in Plucker space

        The sign of the result provides the relative orientation of the two lines (0: intersect, <0: skew ccw, >0: skew cw)
        */
        S dot(const MathPlucker6<S>& p2) const
        {
            return mDirection.dot(p2.mLocation) + mLocation.dot(p2.mDirection);
        }

        /** @brief Return a normalized Plucker point by projection on the unit hypersphere */
        MathPlucker6<S> getNormalized() const
        {
            V_ASSERT(MathArithmetic<S>::isFinite(getDirection()));
            V_ASSERT(MathArithmetic<S>::isFinite(getLocation()));

            MathPlucker6<S>  result(*this);

            S myDiv = sqrt(mDirection.getSquaredNorm() + mLocation.getSquaredNorm());
            V_ASSERT(MathArithmetic<S>::getAbs(myDiv) > MathArithmetic<S>::Tolerance());

            result.mDirection /= myDiv;
            result.mLocation /= myDiv;

            V_ASSERT(MathArithmetic<S>::isFinite(result.getDirection()));
            V_ASSERT(MathArithmetic<S>::isFinite(result.getLocation()));

            return result;
        }

        S getDistance(const MathPlucker6<S>& other) const
        {
            S myDiv = MathVector3_<S>::cross(getDirection(),other.getDirection()).getNorm();

            if (MathArithmetic<S>::getAbs(myDiv) > MathArithmetic<S>::Tolerance())
            {
                return MathArithmetic<S>::getAbs(dot(other)) / myDiv;
            }
            else
            {
                // parallel
                S l1l1 = 1.0 / getDirection().getSquaredNorm();
                S l2 = other.getDirection().getNorm();

                S s = l2 * sqrt(l1l1);
                MathVector3_<S> l = other.getLocation();
                l /= s;
                return MathVector3_<S>::cross(getDirection(),getLocation() - l).getNorm() * l1l1;
            }
        }


        /**@brief Origin point in Plucker space*/
        static MathPlucker6<S>  Zero()
        {
            return MathPlucker6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        }

        friend std::ostream& operator<<(std::ostream& os, const MathPlucker6& v)
        {
            os << v.mDirection << v.mLocation;
            return os;
        }

    private:
        MathVector3_<S> mDirection;   /**< @brief Direction part of the Plucker point*/
        MathVector3_<S> mLocation;    /**< @brief Location part of the Plucker point*/

    };
}
