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

#include <unordered_map>
#include <vector>
#include <math.h>
#include <stack>

#include "math_vector_2.h"
#include "math_vector_3.h"
#include "math_plane_3.h"
#include "math_arithmetic.h"
#include "math_geometry.h"
#include "math_predicates.h"
#include "geometry_convex_polygon.h"
#include "silhouette_mesh_face.h"
#include "geometry_convex_hull.h"
#include "helper_debug_visualisation.h"
#include "helper_statistic_collector.h"
#include "silhouette.h"
#include "visilib_core.h"

namespace visilib
{

    /** @brief Provides the implemention of the silhouette optimization heuristics

    Instead of splitting the polytopes by all the edges of the occluders, the silhouette optimization heuristic detects and discards all the internal edges,
    that do not participate to the visibility between the source polygons. This is performed by extracting the silhouette edges of a triangle mesh occluder with respect to the
    source polygons. The silhouette detection starts from a face of the mesh. All the faces of the triangle mesh are visited successively by using the adjacency information between the faces.
    During the visit, the potential silhouette edges are detected and added to the silhouette. The resulting silhouette is attached to all the visited faces for later use.
    */

    class SilhouetteProcessor
    {
    public:
        SilhouetteProcessor(HelperStatisticCollector* aStatisticCollector);

        /** @brief Attach a debbuger to store debugging information of the silhouette computations */
        void attachVisualisationDebugger(VisualDebugger* aDebugger) { mDebugger = aDebugger; }

        /** @brief Initialize the silhouette processor for the two source polygons */
        void init(const GeometryConvexPolygon& aSource1, const GeometryConvexPolygon& aSource2);

        /** @brief Extract all the silhouettes with respect to the two source polygons

        The computed silhouette are attached to all the visited faces during silhouette extraction.
        */
        void extractSilhouette(size_t geometryId, const std::vector<SilhouetteMeshFace>& faces, bool silhouetteOptimization, std::vector< VisibilitySilhouette*>& silhouettes);

        /** @brief Find the silhouette associated to a given face
        @param face: the index of the face we are looking for a silhouette
        @return: the silhouette if it exists, nullptr otherwise
        */
        VisibilitySilhouette* findSilhouette(SilhouetteMeshFace* face)
        {
            auto iter = mSilhouetteCache.find(face);

            return iter == mSilhouetteCache.end() ? nullptr : iter->second;
        }

    private:

        /** @brief Test if the edge joining two faces is potentially a silhouette edge with respect to the two source polygons */
        bool isPotentialSilhouetteEdge(SilhouetteMeshFace* face0, SilhouetteMeshFace* face1);

        /** @brief  Actual implementation of the Test if the edge joining two faces is potentially a silhouette edge with respect to the two source polygons */
        bool isPotentialSilhouetteEdgeInternal(SilhouetteMeshFace* face0, SilhouetteMeshFace* face1);

        /** @brief Test if the provided face lies between the two source polygons */
        bool isPolygonBetweenSourcePlanes(SilhouetteMeshFace* face);

        /** @brief Compute the convex hull between the two souce polygons.   */
        void initConvexHull();

        /** @brief Add a silhouette to a face*/
        void addSilhouette(SilhouetteMeshFace* face, VisibilitySilhouette* silhouette)
        {
            mSilhouetteCache.insert(std::pair<SilhouetteMeshFace*, VisibilitySilhouette*>(face, silhouette));
        }

        const GeometryConvexPolygon* mSource[2];                                       /**< @brief The two source polygons*/
        std::unordered_map<SilhouetteMeshFace*, VisibilitySilhouette*>mSilhouetteCache;  /**< @brief A silhouette cache associating a silhouette to a face of an occluder*/
        std::unordered_map<SilhouetteMeshFace*, bool> mPolygonBetweenSourcePlanesCache;  /**< @brief A cache storing if a polygon is between the source polygons*/
        std::unordered_map<size_t, bool> mPotentialSilhouetteEdgeCache;                /**< @brief A cache storing if an edge is a potential silhouette edge*/

        GeometryConvexHull* mConvexHull;                                               /**< @brief The convex hull of the source polygons*/
        VisualDebugger* mDebugger;                                                     /**< @brief A debugger collecting debugging information during silhouette extraction*/
        HelperStatisticCollector* mHelperStatisticCollector;                           /**< @brief A statistic collector collecting statistics during silhouette extraction */
    };

    inline SilhouetteProcessor::SilhouetteProcessor(HelperStatisticCollector* aHelperStatisticCollector)
        : mDebugger(nullptr),
        mConvexHull(nullptr),
        mHelperStatisticCollector(aHelperStatisticCollector)
    {
    }

    inline void SilhouetteProcessor::init(const GeometryConvexPolygon& aSource1, const GeometryConvexPolygon& aSource2)
    {
        mSource[0] = &aSource1;
        mSource[1] = &aSource2;
        initConvexHull();
    }

    inline void SilhouetteProcessor::initConvexHull()
    {
        mConvexHull = GeometryConvexHullBuilder::build(mSource[0]->getVertices(), mSource[1]->getVertices());

        V_ASSERT(mConvexHull != nullptr);
    }

    //--------------------------------------------------------------------------------------
    // All silhouette optimization code from Otso Makinen
    // Hybrid Software and Umbra Software
    //--------------------------------------------------------------------------------------

    inline bool SilhouetteProcessor::isPolygonBetweenSourcePlanes(SilhouetteMeshFace * face)
    {
        auto iter = mPolygonBetweenSourcePlanesCache.find(face);
        if (iter != mPolygonBetweenSourcePlanesCache.end())
            return iter->second;

        bool inside = true;

        for (int i = 0; i < 2; i++)
        {
            bool vertexOnFrontSide = MathGeometry::hasVertexOnFrontSide(mSource[i]->getPlane(), face);

            if (!vertexOnFrontSide)
            {
                inside = false;
            }
        }
        mPolygonBetweenSourcePlanesCache.insert(std::pair<SilhouetteMeshFace*, bool>(face, inside));
        return inside;
    }

    static inline size_t getKey(SilhouetteMeshFace * face0, SilhouetteMeshFace * face1)
    {
        return ((size_t)(face0)+(size_t)(face1));
    }

    inline bool SilhouetteProcessor::isPotentialSilhouetteEdge(SilhouetteMeshFace * face0, SilhouetteMeshFace * face1)
    {
        size_t key = getKey(face0, face1);
        auto iter = mPotentialSilhouetteEdgeCache.find(key);
        if (iter != mPotentialSilhouetteEdgeCache.end())
            return iter->second;

        bool result = isPotentialSilhouetteEdgeInternal(face0, face1);
        mPotentialSilhouetteEdgeCache.insert(std::pair<size_t, bool>(key, result));
        return result;
    }

    inline static bool intersect(const GeometryConvexPolygon & polygon, const MathPlane3d & aPlane)
    {
        return MathPredicates::getRelativePosition(polygon.getVertices(), aPlane) == ON_BOUNDARY;
    }

    inline bool SilhouetteProcessor::isPotentialSilhouetteEdgeInternal(SilhouetteMeshFace * face0, SilhouetteMeshFace * face1)
    {
        if (!isPolygonBetweenSourcePlanes(face0) || !isPolygonBetweenSourcePlanes(face1))
        {
            return false;
        }

        MathPlane3d plane0 = convert<MathPlane3d>(MathGeometry::computePlane(face0->getVertex(0), face0->getVertex(1), face0->getVertex(2)));
        MathPlane3d plane1 = convert<MathPlane3d>(MathGeometry::computePlane(face1->getVertex(0), face1->getVertex(1), face1->getVertex(2)));

        for (int s = 0; s < 2; s++)
        {
            if (intersect(*mSource[s], plane0) || intersect(*mSource[s], plane1))
            {
                continue;
            }
            const MathVector3d& myVertex = mSource[s]->getVertex(0);
            double dot0 = plane0.dot(myVertex);
            double dot1 = plane1.dot(myVertex);
            if (MathArithmetic<double>::isSameSign(dot0, dot1))
            {
                return false;
            }
        }

        if (!MathArithmetic<double>::isSameSign(plane0.d, plane1.d))
        {
            return true;
        }

        if (MathGeometry::hasVertexOutsidePlane(plane1, face0))
        {
            return true;
        }

        return false;
    }

    inline void SilhouetteProcessor::extractSilhouette(size_t geometryId, const std::vector<SilhouetteMeshFace> & meshFaces, bool silhouetteOptimization, std::vector<VisibilitySilhouette*> & silhouettes)
    {
        std::vector<bool> processed;
        processed.resize(meshFaces.size());

        for (size_t faceIndex = 0; faceIndex < meshFaces.size(); faceIndex++)
        {
            if (processed[faceIndex])
                continue;

            std::stack<int> stack;
            stack.push((int)faceIndex);

            VisibilitySilhouette* s = new VisibilitySilhouette(meshFaces, geometryId);
            silhouettes.push_back(s);
            
            while (!stack.empty())
            {
                int myIndex = stack.top();                
                stack.pop();
                if (processed[myIndex])
                    continue;

                V_ASSERT(!processed[myIndex]);
                processed[myIndex] = true;
                SilhouetteMeshFace* face = const_cast<SilhouetteMeshFace*>(&meshFaces[myIndex]);
                V_ASSERT(findSilhouette(face) == nullptr);

                if (mConvexHull != nullptr)
                {
                    std::vector<bool> hasNeighbours(face->getVertexCount(), true);

                    bool faceIsInsideHull = false;
                    for (int edgeIndex = 0; edgeIndex < face->getVertexCount(); edgeIndex++)
                    {
                        MathVector2i edge = face->getEdge(edgeIndex);

                        if (MathGeometry::isEdgePotentiallyInsideShaft(mConvexHull->getFaces(), face->getVertex(edge.x), face->getVertex(edge.y), false))
                        {
                            faceIsInsideHull = true;

                            if (!silhouetteOptimization 
                                || face->getNeighbours(edgeIndex) == -1 || isPotentialSilhouetteEdge(face, const_cast<SilhouetteMeshFace*>(&meshFaces[face->getNeighbours(edgeIndex)])))
                            {
                                s->addEdge(face, edgeIndex,0);
                                hasNeighbours[edgeIndex] = false;
                            }
                        }
                        else
                        {
                            hasNeighbours[edgeIndex] = false;
                        }
                    }

                    if (faceIsInsideHull)
                    {
                        addSilhouette(face, s);

                        if (mDebugger != nullptr)
                        {
                            mDebugger->get()->addRemovedTriangle(face->getVertex(0), face->getVertex(1), face->getVertex(2));
                        }
                        
                        s->addFace(*face);

                        mHelperStatisticCollector->inc(OCCLUDER_TRIANGLE_COUNT);

                        if (silhouetteOptimization)
                        {
                            for (int i = 0; i < face->getVertexCount(); i++)
                            {
                                if (hasNeighbours[i]&&!processed[face->getNeighbours(i)])
                                {
                                    stack.push(face->getNeighbours(i));
                                }
                            }
                        }

                    }
                }
            }
         }
    }
}
