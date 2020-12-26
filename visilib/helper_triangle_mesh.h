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

#include <vector>
#include <iostream>

#include "geometry_aabbox.h"
#include "math_vector_3.h"
#include "math_arithmetic.h"

namespace visilib
{
    /** @brief Store an indexed triangle mesh in 3D space */
    class HelperTriangleMesh
    {
    public:
        HelperTriangleMesh(const std::vector<float>& vertices, const std::vector<int>& triangles);
        HelperTriangleMesh(const std::vector<MathVector3f>& vertices, const std::vector<int>& triangles);
        ~HelperTriangleMesh();

        /** @brief Compute additional information (bounding box and normals) */
        void prepare();

        size_t getVertexCount() const
        {
            return mVertices.size();
        }

        size_t getIndiceCount() const
        {
            return mIndices.size();
        }

        const std::vector<MathVector3f>& getVertices() const
        {
            return mVertices;
        }

        std::vector<MathVector3f>& getVertices()
        {
            return mVertices;
        }

        const std::vector<int>& getIndices() const
        {
            return mIndices;
        }

        std::vector<int>& getIndices()
        {
            return mIndices;
        }

        const std::vector<MathVector3f>& getNormals() const
        {
            return mNormals;
        }

        /** @brief Remove the degenerated triangles of the mesh
        
        A triangle is degenerated if one of its edge is smaller than the provided threshold
        @param threshold: the threshold used to test the edges length
        */
        void removeDegeneratedTriangles(float threshold);

    private:

        /** @brief Return the number of triangles of the triangle mesh */
        size_t	getTriangleCount(){ return (size_t)mIndices.size() / 3; };

        /** @brief Compute and store the vertex normals of the triangle mesh */
        void computeNormal();

        void swapLastIndice(size_t i);

        std::vector<int>   mIndices;         /** < @brief The index table of the mesh */
        GeometryAABB mBox;                   /** < @brief The axis aligned bounding box of the mesh*/
        std::vector<MathVector3f> mVertices; /** < @brief The list of vertices of the mesh in 3D space */
        std::vector<MathVector3f> mNormals;  /** < @brief The list of vertex normal in 3D space */
    };

    inline HelperTriangleMesh::HelperTriangleMesh(const std::vector<float>& vertices, const std::vector<int>& triangles)
    {
        mIndices = triangles;
        mVertices.resize(vertices.size() / 3);

        for (size_t i = 0; i < mVertices.size(); i++)
        {
            mVertices[i] = MathVector3f(vertices[i * 3], vertices[i * 3 + 1], vertices[i * 3 + 2]);
        }
    }

    inline HelperTriangleMesh::HelperTriangleMesh(const std::vector<MathVector3f>& vertices, const std::vector<int>& triangles)
    {
        mIndices = triangles;
        mVertices = vertices;
    }

    inline void HelperTriangleMesh::prepare()
    {
        mNormals.resize(mVertices.size());
        computeNormal();
    }

    inline HelperTriangleMesh::~HelperTriangleMesh()
    {
    }

    inline void HelperTriangleMesh::computeNormal()
    {
        if (mIndices.size() <= 4)
        {
            MathVector3f d1 = mVertices[1] - mVertices[0];
            MathVector3f d2 = mVertices[2] - mVertices[0];

            mNormals[0] = MathVector3f::cross(d1, d2);
            mNormals[0].normalize();
            mNormals[1] = mNormals[0];
            mNormals[2] = mNormals[0];
            return;
        }

        for (size_t i = 0; i < mVertices.size(); i++)
        {
            mNormals[i] = MathVector3f::Zero();
        }
        for (size_t i = 0; i < mIndices.size() - 4; i = i + 3)
        {
            MathVector3f d1 = mVertices[mIndices[i + 1]] - mVertices[mIndices[i]];
            MathVector3f d2 = mVertices[mIndices[i + 2]] - mVertices[mIndices[i]];

            MathVector3f n = MathVector3f::cross(d1, d2);
            n.normalize();
            mNormals[mIndices[i]] += n;
            mNormals[mIndices[i + 1]] += n;
            mNormals[mIndices[i + 2]] += n;
        }

        for (size_t i = 0; i < (size_t)mVertices.size(); i++)
        {
            mNormals[i].normalize();
        }
    }

    inline void HelperTriangleMesh::swapLastIndice(size_t i)
    {
        if ((size_t)mIndices.size() - 3>i + 2)
        {
            mIndices[i] = mIndices[mIndices.size() - 3];
            mIndices[i + 1] = mIndices[mIndices.size() - 2];
            mIndices[i + 2] = mIndices[mIndices.size() - 1];
        }
        mIndices.resize(mIndices.size() - 3);
    }

    inline void HelperTriangleMesh::removeDegeneratedTriangles(float threshold)
    {
        MathVector3f myNormal;
        MathVector3f d1, d2;
        if (mIndices.size() == 0)
        {
            std::cout << "Warning: Mesh without triangle" << std::endl;
            return;
        }

        if (mVertices.size() == 0)
        {
            std::cout << "Warning: Mesh without vertex" << std::endl;
            return;
        }

        if (mIndices.size() <= 4)
        {
            myNormal = MathVector3f::cross(mVertices[1] - mVertices[0], mVertices[2] - mVertices[0]);
            if (myNormal.isZero(threshold))mIndices.clear();
        }
        else
        {
            for (size_t i = 0; i < mIndices.size() - 3; i += 3)
            {
                V_ASSERT(mIndices[i + 0] < mVertices.size());
                V_ASSERT(mIndices[i + 1] < mVertices.size());
                V_ASSERT(mIndices[i + 2] < mVertices.size());
                myNormal = MathVector3f::cross(mVertices[mIndices[i + 1]] - mVertices[mIndices[i]], mVertices[mIndices[i + 2]] - mVertices[mIndices[i]]);
                if (myNormal.isZero(threshold))
                {
                    swapLastIndice(i);
                    i -= 3;
                }
            }
        }
    }
}
