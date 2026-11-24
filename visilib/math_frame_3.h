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
#include "math_vector_2.h"
#include "math_vector_3.h"
#include "math_plane_3.h"

namespace visilib
{
    template<class S>
    class MathFrame3_;

    typedef MathFrame3_<float> MathFrame3f;
    typedef MathFrame3_<double> MathFrame3d;

    /** @brief An orthonormal frame of a plane: an origin, two unit directions spanning the plane, and the unit normal.

    The directions satisfy mAxisU x mAxisV = mNormal. A point of the plane has the coordinates (u, v) with
    p = mOrigin + u mAxisU + v mAxisV; the frame converts between the two representations, and gives the
    plane equation n . x + d = 0 it lies on.
        @param S the coordinates representation (float, double, exact)*/
    template<class S>
    class MathFrame3_
    {
    public:
        MathFrame3_()
            : mOrigin(MathVector3_<S>::Zero()), mAxisU(MathVector3_<S>::X()), mAxisV(MathVector3_<S>::Y()), mNormal(MathVector3_<S>::Z())
        {
        }

        /** @brief Build the frame of the plane through anOrigin with the normal aNormal, aFirstAxis giving the direction of u.

        aFirstAxis is projected onto the plane: it needs only not to be parallel to the normal. The normal and the
        axis need not be unit vectors.*/
        MathFrame3_(const MathVector3_<S>& anOrigin, const MathVector3_<S>& aNormal, const MathVector3_<S>& aFirstAxis)
            : mOrigin(anOrigin), mNormal(aNormal)
        {
            mNormal.normalize();
            mAxisU = aFirstAxis - mNormal * mNormal.dot(aFirstAxis);
            mAxisU.normalize();
            mAxisV = MathVector3_<S>::cross(mNormal, mAxisU);
        }

        /** @brief The coordinates (u, v) of the projection of aPoint onto the plane.*/
        MathVector2_<S> toLocal(const MathVector3_<S>& aPoint) const
        {
            const MathVector3_<S> d = aPoint - mOrigin;
            return MathVector2_<S>(d.dot(mAxisU), d.dot(mAxisV));
        }

        /** @brief The point of the plane with the coordinates (u, v).*/
        MathVector3_<S> toWorld(const MathVector2_<S>& aPoint) const
        {
            return mOrigin + mAxisU * aPoint.x + mAxisV * aPoint.y;
        }

        /** @brief The signed distance of aPoint to the plane, positive on the side of the normal.*/
        S getSignedDistance(const MathVector3_<S>& aPoint) const
        {
            return (aPoint - mOrigin).dot(mNormal);
        }

        /** @brief The plane n . x + d = 0 of the frame.*/
        MathPlane3_<S> getPlane() const
        {
            return MathPlane3_<S>(mNormal, -mNormal.dot(mOrigin));
        }

        const MathVector3_<S>& getOrigin() const { return mOrigin; }
        const MathVector3_<S>& getAxisU() const { return mAxisU; }
        const MathVector3_<S>& getAxisV() const { return mAxisV; }
        const MathVector3_<S>& getNormal() const { return mNormal; }

    private:
        MathVector3_<S> mOrigin;  /**< @brief The point of coordinates (0, 0)*/
        MathVector3_<S> mAxisU;   /**< @brief The unit direction of the u coordinate*/
        MathVector3_<S> mAxisV;   /**< @brief The unit direction of the v coordinate, mNormal x mAxisU*/
        MathVector3_<S> mNormal;  /**< @brief The unit normal of the plane*/
    };
}
