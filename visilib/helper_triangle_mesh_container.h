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

#include "visilib.h"

#include <vector>
#include <unordered_map>
#include "helper_triangle_mesh.h"
#include "geometry_ray.h"
#include "math_geometry.h"
#include "geometry_aabbox.h"

namespace visilib
{
    class HelperTriangleMesh;


    /** @brief Represents a scene containing all the geometry represented as GeometeryMesh
    */

    class HelperTriangleMeshContainer
    {
    public:
        HelperTriangleMeshContainer();

        virtual ~HelperTriangleMeshContainer();


        /** @brief Return the list of ids of the triangle meshes stored in the scene*/
        const std::vector<size_t>& getGeometryIds() const
        {
            return mGeometryIds;
        }

        size_t getGeometryId(size_t index) const
        {
            return mGeometryIds[index];
        }

        size_t getGeometryCount() const
        {
            return mGeometryIds.size();
        }

        /** @brief Add a mesh to the scene*/
        virtual size_t add(HelperTriangleMesh* mesh);

        /** @brief Return triangle mesh data information of a geometry stored in the scene*/
        GeometryDiscreteMeshDescription* createTriangleMeshDescription(size_t geometryId) const;

        /** @brief Remove the degnerated triangles*/
        virtual void removeDegeneratedTriangles(float threshold);

        /** @brief Return the list of tirangle meshes of the scene */
        std::vector<HelperTriangleMesh*>& getMeshArray()
        {
            return mMeshArray;
        }

    protected:
     
        std::vector<HelperTriangleMesh*> mMeshArray;  /**< @brief The list of triangle meshes representing the occluders set*/
        std::vector<size_t> mGeometryIds;       /**< @brief The ids of the triangle meshes*/
    };

    inline HelperTriangleMeshContainer::HelperTriangleMeshContainer()
    {
    }

    inline HelperTriangleMeshContainer::~HelperTriangleMeshContainer()
    {
        for (size_t i = 0; i < mMeshArray.size(); i++)
        {
            delete mMeshArray[i];
        }
    }

    inline size_t HelperTriangleMeshContainer::add(HelperTriangleMesh* mesh)
    {
        size_t id = mMeshArray.size();
        mGeometryIds.push_back(id);
        mMeshArray.push_back(mesh);
        mesh->prepare();
        return id;
    }

    inline void HelperTriangleMeshContainer::removeDegeneratedTriangles(float threshold)
    {
        for (auto iter = mMeshArray.begin(); iter != mMeshArray.end(); iter++)
        {
            (*iter)->removeDegeneratedTriangles(threshold);
        }
    }

    inline GeometryDiscreteMeshDescription* HelperTriangleMeshContainer::createTriangleMeshDescription(size_t geometryId) const
    {
        HelperTriangleMesh* mesh = mMeshArray[geometryId];
        GeometryTriangleMeshDescription* info = new GeometryTriangleMeshDescription;
        info->vertexCount = mesh->getVertexCount();
        info->faceCount= mesh->getIndiceCount() / 3;
        info->indexArray = &mesh->getIndices()[0];
        info->vertexArray = (float*)& mesh->getVertices()[0];

        info->normalArray = mesh->getNormals().empty() ? nullptr : (float*)& mesh->getNormals()[0];
        return info;
    }
}
