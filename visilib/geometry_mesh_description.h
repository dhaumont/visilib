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

namespace visilib
{

    /** @brief Description of a discrete geometry (ie a mesh).*/

    struct GeometryDiscreteMeshDescription
    {
        GeometryDiscreteMeshDescription();

        virtual ~GeometryDiscreteMeshDescription()
        {
        }

        virtual std::vector<int> getIndices(size_t aFace) const = 0;
        virtual size_t getIndexCount() const = 0;

        size_t vertexCount;         /**< @brief Number of vertices*/
        size_t faceCount;           /**< @brief Number of faces*/
        const float* vertexArray;   /**< @brief Pointer to the vertex table*/
        const float* normalArray;   /**< @brief Pointer to the normal table (optional)*/

    };

    /** @brief Representation of a triangle mesh.*/

    struct GeometryTriangleMeshDescription : public GeometryDiscreteMeshDescription
    {
        GeometryTriangleMeshDescription();

        virtual std::vector<int> getIndices(size_t aFace) const override
        {
            return std::vector<int>{indexArray[aFace * 3], indexArray[aFace * 3 + 1], indexArray[aFace * 3 + 2] };
        }

        virtual size_t getIndexCount() const override
        {
            return faceCount * 3;
        }

        const int* indexArray;      /**< @brief Pointer to the index table*/
    };

inline GeometryDiscreteMeshDescription::GeometryDiscreteMeshDescription()
    : vertexCount(0),
    faceCount(0),
    vertexArray(nullptr),
    normalArray(nullptr)
{
}

inline GeometryTriangleMeshDescription::GeometryTriangleMeshDescription()
    : GeometryDiscreteMeshDescription(),
    indexArray(nullptr)
{
}
}
