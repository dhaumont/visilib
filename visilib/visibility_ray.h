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

namespace visilib
{
    /** @brief Represents a 3D ray, and contains the list of intersected geometries.*/

    struct VisibilityRay
    {
        float org[3];          /**< @brief Ray origin*/
        float dir[3];          /**< @brief Ray direction*/
        float tnear;           /**< @brief Start of ray segment*/
        float tfar;            /**< @brief End of ray segment (set to hit distance)*/

        void addIntersection(size_t geometryId, size_t primitiveId, double t)
        {
            mGeometryIds.push_back(geometryId);
            mPrimitiveIds.push_back(primitiveId);
            mTs.push_back(t);
        }

        std::vector<size_t> mGeometryIds;     /**< @brief Intersected geometry mesh*/
        std::vector<size_t> mPrimitiveIds;    /**< @brief Intersected primitve of the mesh*/
        std::vector<double> mTs;    /**< @brief Intersected primitve of the mesh*/

        VisibilityRay() {
            tnear = 0.0f;
            tfar = 1.0f;

            dir[0] = dir[1] = dir[2] = 0.0f;
            org[0] = org[1] = org[2] = 0.0f;
        }
    };
}
