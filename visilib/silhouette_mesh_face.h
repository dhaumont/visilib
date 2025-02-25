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
#include "visilib.h"
#include "math_vector_2.h"
#include "math_vector_3.h"
namespace visilib
{
    /** < @brief Storage of the connectivity between mesh faces of the occluder geometry. Each occluder is stored as a set of SilhouetteMeshFace. 
    Each SilhouetteMeshFace contains the definition of the vertices of the face as well as a pointer to the neigbour faces.
    */

    class SilhouetteMeshFace
    {
    public:
        SilhouetteMeshFace(size_t n = 3)
        {
            mVertices.resize(n,nullptr);
            mNeighbours.resize(n,-1);
 
            mFaceIndex = 0; 
        }

        SilhouetteMeshFace(const SilhouetteMeshFace& other)
        {
            mVertices = other.mVertices;
            mNeighbours = other.mNeighbours;

            mFaceIndex = other.mFaceIndex;
        }

        ~SilhouetteMeshFace()
        {
        }

        /** @brief Initialize a face from a triangle mesh*/
        void setGeometry(GeometryDiscreteMeshDescription* aMesh, size_t aFace)
        {
            const MathVector3f* myVertices = (MathVector3f*)aMesh->vertexArray;
            std::vector<int> myIndices = aMesh->getIndices(aFace);
 
            for (size_t i = 0; i < myIndices.size(); i++)
            {
                mVertices[i] = (MathVector3f*)& myVertices[myIndices[i]];
            }
            mFaceIndex = aFace;
        }

        /** @brief Store the neigbour of a face, sharing the edge i */

        void setNeighbour(size_t i, int neighbour)
        {
            mNeighbours[i] = neighbour;
        }

        /** @brief Return the neigbour of the face sharing the edge i */

        int getNeighbours(size_t i) const
        {
            return mNeighbours[i];
        }

        /** @brief Return a vertex of the face */

        const MathVector3f& getVertex(size_t i) const
        {
            return *mVertices[i];
        }

        size_t getVertexCount() const
        {
            return mVertices.size();
        }

        size_t getFaceIndex() const
        {
            return mFaceIndex;
        }

        /** @brief Return a 2D segment joining the two vertices of an edge  */

        MathVector2i getEdge(size_t i) const
        {
            return MathVector2i(i, (i + 1) % getVertexCount());
        }

    private:
        std::vector<MathVector3f*> mVertices;
        std::vector<int> mNeighbours;
        size_t mFaceIndex;
    };
}
