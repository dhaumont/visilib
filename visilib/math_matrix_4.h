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

#include <string.h> //for memset
#include <math.h>

#include "math_vector_3.h"

namespace visilib
{
    template<class S>
    class MathMatrix_;
    typedef MathMatrix_<float> MathMatrixf;
    typedef MathMatrix_<double> MathMatrixd;

    /** @brief Represents a 4x4 transformation matrix in 3D space.
       @param S the coordinates representation (float, double, exact)
    */
    template<class S>
    class MathMatrix_
    {
    public:
        MathMatrix_()
        {
            memset(&f[0], 0, 16 * sizeof(S));
        }

        MathVector3_<S>  multiply(const MathVector3_<S>& aVector) const
        {
            return MathVector3_<S>(
                v[0][0] * aVector.x + v[1][0] * aVector.y + v[2][0] * aVector.z + v[3][0],
                v[0][1] * aVector.x + v[1][1] * aVector.y + v[2][1] * aVector.z + v[3][1],
                v[0][2] * aVector.x + v[1][2] * aVector.y + v[2][2] * aVector.z + v[3][2]);
        }

        MathPlane3_<S> multiply(const MathPlane3_<S> & aVector) const
        {
            return MathPlane3_<S>(
                v[0][0] * aVector.mNormal.x + v[1][0] * aVector.mNormal.y + v[2][0] * aVector.mNormal.z + v[3][0] * aVector.d,
                v[0][1] * aVector.mNormal.x + v[1][1] * aVector.mNormal.y + v[2][1] * aVector.mNormal.z + v[3][1] * aVector.d,
                v[0][2] * aVector.mNormal.x + v[1][2] * aVector.mNormal.y + v[2][2] * aVector.mNormal.z + v[3][2] * aVector.d,
                v[0][3] * aVector.mNormal.x + v[1][3] * aVector.mNormal.y + v[2][3] * aVector.mNormal.z + v[3][3] * aVector.d);
        }

        void setRotateZ(S phi)
        {
            memset(&f[0], 0, 16 * sizeof(S));
            v[2][2] = 1; v[3][3] = 1;

            //rotate z
            v[0][0] = cos(phi);
            v[0][1] = sin(phi);

            v[1][0] = -sin(phi);
            v[1][1] = cos(phi);
        }

        void  setRotateY(S alpha)
        {
            memset(&f[0], 0, 16 * sizeof(S));
            v[1][1] = 1; v[3][3] = 1;

            //rotate y
            v[0][0] = cos(alpha);
            v[0][2] = -sin(alpha);
            v[2][0] = sin(alpha);
            v[2][2] = cos(alpha);
        }

        void  setRotateX(S theta)
        {
            memset(&f[0], 0, 16 * sizeof(S));
            v[0][0] = 1; v[3][3] = 1;

            //rotate x
            v[1][1] = cos(theta);
            v[1][2] = sin(theta);
            v[2][1] = -sin(theta);
            v[2][2] = cos(theta);
        }

        void setTranslation(const MathVector3_<S> & aCenter)
        {
            v[3][0] = aCenter.x;
            v[3][1] = aCenter.y;
            v[3][2] = aCenter.z;
        }

    private:
        union
        {
            S v[4][4];
            S f[16];
        };
    };
}
