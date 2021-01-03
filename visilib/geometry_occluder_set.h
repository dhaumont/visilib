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

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <stack>
#include <queue>
#include "silhouette_mesh_face.h"
#include "geometry_ray.h"
#include "math_geometry.h"

namespace visilib
{
    struct GeometryDiscreteMeshDescription;
    class SilhouetteMeshFace;
   
    /** @brief Stores the occluders against which visibility is tested. The occluders are stored under the form of a connected set of faces, that are used for efficient silhouette detection.
    Connectivity information of the occluders is computed in a lazy way, only when required.*/

    class GeometryOccluderSet
    {
    public:
        void addOccluder(GeometryDiscreteMeshDescription* info);

        /** @brief Prepare the scene before ray tracing */
        void prepare();

        size_t getOccluderCount()const
        {
            return mOccluders.size();
        }

        /** @brief Return the list of connected faces of a mesh

        @param scene: the scene containing the triangle mesh
        @param geometryId: the index of the triangle mesh
        @return : the list of faces of the mesh, each facces contains a links to its neighbour faces
        */

        std::vector<SilhouetteMeshFace>* getOccluderConnectedFaces(size_t geometryId);

        /** @brief Restore the faces of the scene to their initial state, discarding any clipping if it has been performed*/

        void restoreOccluderConnectedFaces();

        /** @brief Store the list of faces of a mesh

        @param mesh: the triangle mesh
        @param: the list of faces
        */
        void setOccluderConnectedFaces(GeometryDiscreteMeshDescription* mesh, std::vector<SilhouetteMeshFace>& aFaces);

        ~GeometryOccluderSet()
        {
            for (auto iter : mConnectedFacesCache)
            {
                delete iter;
            }

            for (auto iter : mOccluders)
            {
                delete iter;
            }

        }
    private:

        /** @brief Provides a unique key for each edge of the triangle mesh. The edge is ordered*/

        size_t getKey(size_t v0, size_t v1, size_t myFaceNumber)
        {
            return v0 < v1 ? v0 * myFaceNumber + v1 : v1 * myFaceNumber + v0;
        }

        /** @brief Compute the list of faces of a triangle mesh, containing the adjacency information

        The adjacency is computed in a single pass over the triangle mesh, by constructing on the fly the edge connectivity information

        @param mesh: the triangle mesh
        @param: the list of faces
        */
        void extractConnectedMeshFaces(GeometryDiscreteMeshDescription* mesh, std::vector<SilhouetteMeshFace>& aFaces);


        /**@brief The list of faces of the triangle meshes.

        The faces are stored in a map indexed by the id of the mesh
        */
        std::vector<std::vector<SilhouetteMeshFace>*> mConnectedFacesCache;
        std::unordered_map<size_t, size_t> mLastHit;
        std::vector<GeometryDiscreteMeshDescription*> mOccluders;
        std::vector<GeometryAABB> mBoundingBoxes;
    };

    inline std::vector<SilhouetteMeshFace>* GeometryOccluderSet::getOccluderConnectedFaces(size_t geometryId)
    {
        std::vector<SilhouetteMeshFace>* myFaces;
      
        if (mConnectedFacesCache[geometryId] == nullptr)
        {  
            GeometryDiscreteMeshDescription* mesh = mOccluders[geometryId];
     
            myFaces = new std::vector<SilhouetteMeshFace>();
            extractConnectedMeshFaces(mesh, *myFaces);
            mConnectedFacesCache[geometryId] = myFaces;
        }
        return mConnectedFacesCache[geometryId];
    }

    inline void GeometryOccluderSet::restoreOccluderConnectedFaces()
    {
        //remove the clipping performed if any
        for (size_t geometryId = 0; geometryId < mConnectedFacesCache.size(); geometryId++)
        {
            GeometryDiscreteMeshDescription* mesh = mOccluders[geometryId];
          
            setOccluderConnectedFaces(mesh, *mConnectedFacesCache[geometryId]);
        }
    }

    inline void GeometryOccluderSet::extractConnectedMeshFaces(GeometryDiscreteMeshDescription* mesh, std::vector<SilhouetteMeshFace>& aFaces)
    {
        size_t myFaceNumber = mesh->faceCount;

        aFaces.clear();
        aFaces.resize(myFaceNumber);

        for (size_t i = 0; i < myFaceNumber; i++)
        {
            aFaces[i].setGeometry(mesh, i);
        }
        std::unordered_map<size_t, MathVector2i> myEdgesMap;

        for (size_t myTriangle = 0; myTriangle < myFaceNumber; myTriangle++)
        {
            const std::vector<int> myIndices = mesh->getIndices(myTriangle);

            for (size_t myEdge = 0; myEdge < 3; myEdge++)
            {
                size_t myBegin = myIndices[myEdge];
                size_t myEnd = myIndices[(myEdge + 1) % myIndices.size()];

                const size_t key = getKey(myBegin, myEnd, mesh->faceCount*3);
                auto iter = myEdgesMap.find(key);
                if (iter == myEdgesMap.end())
                {
                    aFaces[myTriangle].setNeighbour(myEdge, -1);
                    myEdgesMap.insert(std::pair<size_t, MathVector2i>(key, MathVector2i(myTriangle, myEdge)));
                }
                else
                {
                    size_t myNeighborTriangle = iter->second.x;
                    size_t myNeighborEdge = iter->second.y;

                    myEdgesMap.erase(iter);

                    aFaces[myTriangle].setNeighbour(myEdge, (int)myNeighborTriangle);
                    aFaces[myNeighborTriangle].setNeighbour(myNeighborEdge, (int)myTriangle);
                }
            }
        }
    }

    inline void GeometryOccluderSet::setOccluderConnectedFaces(GeometryDiscreteMeshDescription* mesh, std::vector<SilhouetteMeshFace> & aFaces)
    {
        size_t myFaceNumber = mesh->faceCount;

        for (size_t i = 0; i < myFaceNumber; i++)
        {
            aFaces[i].setGeometry(mesh, i);
        }
    }

    inline void GeometryOccluderSet::addOccluder(GeometryDiscreteMeshDescription* info)
    {
        mOccluders.push_back(info);
        mConnectedFacesCache.push_back(nullptr);
    }

    /** @brief Prepare the scene before ray tracing */
    inline void GeometryOccluderSet::prepare()
    {
        for (size_t i = 0; i < mOccluders.size(); i++)
        {
            GeometryDiscreteMeshDescription* myTriangleMesh = mOccluders[i];
            MathVector3f myMin, myMax;
            MathArithmetic<float>::getMinMax(myTriangleMesh->vertexArray, myTriangleMesh->vertexCount, myMin, myMax);
            mBoundingBoxes.push_back(GeometryAABB(myMin, myMax));
        }
    }
}
