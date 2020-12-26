/*
Visilib, an open source library for exact visibility computation.
Copyright(C) 2017 by Denis Haumont

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

#include "math_vector_3.h"
#include <iostream>

namespace visilib
{
    template<class S>
    class MathPlane3_;
    typedef MathPlane3_<float> MathPlane3f;
    typedef MathPlane3_<double> MathPlane3d;


    /** @brief Represents a plane in 3D space.
       @param S the coordinates representation (float, double, exact)
    */
    template<class S>
    class MathPlane3_
    {
    public:

        MathPlane3_() : mNormal(), d(0) {	}

        MathPlane3_(S aX, S aY, S aZ, S aD) : mNormal(aX, aY, aZ), d(aD) {}
        MathPlane3_(const MathVector3_<S>& aNormal, S aD) : mNormal(aNormal), d(aD) {}

        void set(const MathVector3_<S>& aNormal, S aW)
        {
            mNormal = aNormal;
            d = aW;
        }

        MathPlane3_<S>& operator += (const MathPlane3_<S>& v) { mNormal += v.mNormal; d += v.d;	return *this; }
        MathPlane3_<S>& operator -= (const MathPlane3_<S>& v) { mNormal -= v.mNormal; d -= v.d;	return *this; }

        static MathPlane3_<S> computeTrianglePlane(const MathVector3_<S>& v0, const MathVector3_<S>& v1, const MathVector3_<S>& v2, S epsilon)
        {
            MathVector3_<S> myNormal = MathVector3_<S>::cross(v2 - v0, v1 - v0);

            S r = myNormal.normalize();

            if (MathArithmetic<S>::getAbs(r) > epsilon)
            {
                return MathPlane3d(myNormal, -myNormal.dot(v0));
            }
            V_ASSERT(0);
            return MathPlane3d::Zero();
        }

        static bool computeTrianglePlane(MathPlane3_<S> & aPlane, MathVector3_<S> * aTriangeVertices, S epsilon)
        {
            aPlane = computeTrianglePlane(aTriangeVertices[0], aTriangeVertices[1], aTriangeVertices[2], epsilon);

            return !aPlane.mNormal.isZero(epsilon);
        }

        S dot(const MathVector3_<S> & aPoint)const
        {
            return mNormal.dot(aPoint) + d;
        }

        static MathPlane3_ Zero() { return MathPlane3_((S)0.0, (S)0.0, (S)0.0, (S)0.0); }

        bool isZero(S anEps)const
        {
            return (mNormal.isZero(anEps) && d <= anEps);
        }

        const MathVector3_<S>& getNormal()const
        {
            return mNormal;
        }
    public:
        MathVector3_<S> mNormal;
        S d;

        friend std::ostream& operator<<(std::ostream & os, const MathPlane3d & plane);
    };

    inline std::ostream& operator<<(std::ostream & os, const MathPlane3d & plane)
    {
        os << "(" << plane.mNormal << "," << plane.d << ")";
        return os;
    }

    template<>
    inline MathPlane3d convert(const MathPlane3f& f)
    {
        return MathPlane3d((double)f.getNormal().x, (double)f.getNormal().y, (double)f.getNormal().z, (double) f.d);
    }
}
