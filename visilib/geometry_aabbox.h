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

#include "math_vector_3.h"

namespace visilib
{
    /** @brief
    Store an axis aligned bounding box in 3D space
    */

    class GeometryAABB
    {
    public:
        GeometryAABB()
            : mMin(0, 0, 0), mMax(0, 0, 0)
        {
        }

        GeometryAABB(const MathVector3f& aMin, const MathVector3f& aMax)
            :mMin(aMin),
            mMax(aMax)
        {
        }

        GeometryAABB(const GeometryAABB& aBox)
            :mMin(aBox.mMin),
            mMax(aBox.mMax)
        {
        }

        const MathVector3f& getMin() const
        {
            return mMin;
        }

        const MathVector3f getMax() const
        {
            return mMax;
        }

        void init(const MathVector3f& aMin, const MathVector3f& aMax)
        {
            mMin = aMin; mMax = aMax;
        }

       
    private:
        MathVector3f mMin;
        MathVector3f mMax;
    };
}