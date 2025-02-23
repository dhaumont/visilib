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

#include <unordered_set>
#include <fstream>
#include "math_plucker_6.h"

namespace visilib
{
    class HelperStatisticCollector;
    template <class P, class S>
    class VisibilityExactQuery_;

    class GeometryConvexPolygon;
    template <class P, class S>
    class PluckerPolytope;
    class SilhouetteMeshFace;
    template <class S>
    class PluckerFacet;
    /** @brief
   Store a silhouette edge, containing the mesh face and the support hyperplane 
   */

    template<class P>
    struct SilhouetteEdge
    {
        SilhouetteEdge(SilhouetteMeshFace* aFace, size_t anEdgeIndex, PluckerFacet<P>* aFacet = nullptr)
            : mEdgeIndex(anEdgeIndex),
            mFace(aFace),
            mPluckerFacet(aFacet),
            mScore(-1.0),
            mIsActive(true)
        {
        }

        PluckerFacet<P>* mPluckerFacet;
        SilhouetteMeshFace* mFace;
        size_t mEdgeIndex;
        double mScore;
        bool mIsActive;
    };


    /** @brief Store a silhouette, representing an occluder as seen from the sources. The silhouettes are stored under the form of a list of silhouette edges. A s
    They are extracted by the silhouette processor class. Each face of the occluder mesh belongs to a single silhouette at any time. 
    In turn, each face of the occluder mesh stores a reference to the silhouette it belongs to.
    */

    template<class P>
    class Silhouette
    {
        /**@brief A hash associated to a silhouette edge for storage in the silhouette */
        struct Hash
        {
        public:

            /**@brief Compute the haskey of a given silhouette edge*/
            std::size_t operator()(std::pair<SilhouetteMeshFace*, size_t> pair) const
            {
                size_t edgeHash = std::hash<size_t>()((size_t)(pair.first));
                size_t faceHash = std::hash<size_t>()(pair.second);

                return faceHash ^ edgeHash;
            }
        };

    public:

        Silhouette(const std::vector<SilhouetteMeshFace>& aFaces, const size_t aGeometryId)
            : mMeshFaces(aFaces),
            mGeometryId(aGeometryId),
            mAvailableEdgeCount(0)
        {
           
        }

        ~Silhouette()
        {
        }

        /**@brief Return the number of edges of the silhouette*/
        int getAvailableEdgeCount() const
        {
            return mAvailableEdgeCount;
        }

        /**@brief Return if the given edge of a given faces belongs to the silhouette*/
   //     bool hasEdge(ConnectedMeshElement* face, size_t edgeIndex) const;

        /**@brief Get the edges of the silhouette*/
        const std::vector<SilhouetteEdge<P>>& getEdges() const { return mEdges; }

        std::vector<SilhouetteEdge<P>>& getEdges() { return mEdges; }

        /**@brief Add an edge to the silhouette*/
        void addEdge(SilhouetteMeshFace* aFace, size_t anEdgeIndex, size_t anHyperplaneIndex);

        void addFace(const SilhouetteMeshFace& face)
        {
            mSilhouetteFaces.push_back(face.getFaceIndex());
        }

        SilhouetteEdge<P>& getEdge(size_t edgeIndex)
        {
            return mEdges[edgeIndex];
        }
        
        const std::vector<size_t>& getEdgesProcessed() { return mEdgesProcessed; }

        void setEdgeActive(size_t edgeIndex, bool isActive)
        {
            SilhouetteEdge<P>& edge = getEdge(edgeIndex);
            V_ASSERT(edge.mIsActive != isActive);

            edge.mIsActive = isActive;
            if (isActive)
                mAvailableEdgeCount++;
            else
                mAvailableEdgeCount--;
            V_ASSERT(mAvailableEdgeCount >= 0 && mAvailableEdgeCount <= mEdges.size());
        }
        void pushEdgeProcessed(size_t edgeIndex)
        {
            mEdgesProcessed.push_back(edgeIndex);
        }

        void popEdgeProcessed(size_t edgeIndex)
        {
            V_ASSERT(mEdgesProcessed.back() == edgeIndex);
            mEdgesProcessed.pop_back();
        }

        const std::vector<size_t>& getSilhouetteFaces() const
        {
            return mSilhouetteFaces;
        }

        const std::vector<SilhouetteMeshFace>& getMeshFaces() const
        {
            return mMeshFaces;
        }

        size_t getGeometryId() const
        {
            return mGeometryId;
        }
    private:
        std::vector<SilhouetteEdge<P>> mEdges;   /**<@brief The list of edges of the silhouette, encoded as pointer to a face and an edge index in the face*/
        std::vector<size_t> mEdgesProcessed;
        std::vector<size_t> mSilhouetteFaces;
        const std::vector<SilhouetteMeshFace>& mMeshFaces;
        const size_t mGeometryId;
        int mAvailableEdgeCount;
    };


    template<class P>
    inline void Silhouette<P>::addEdge(SilhouetteMeshFace* aFace, size_t anEdgeIndex, size_t anHyperplaneIndex)
    {
        mEdges.push_back(SilhouetteEdge<P>(aFace, anEdgeIndex,anHyperplaneIndex));
        mAvailableEdgeCount++;
    }

}

