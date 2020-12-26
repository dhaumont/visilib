/*
Visilib, an open source library for exact visibility computation.
Copyright(C) 2020 by Denis Haumont

This file is part of Visilib.

Visilib is free software : you can redistribute it and / or modify
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

#include <sys/types.h>
#include <vector>
namespace visilib
{
    class VisualDebugger;
    class HelperDebugVisualisation;
    class GeometryOccluderSet;

    /** @brief Represents a 3D ray.*/

    struct Ray
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
    };


    /** @brief Represents a triangle mesh.*/

    struct DiscreteGeometryDescription
    {
        DiscreteGeometryDescription();
        
        virtual ~DiscreteGeometryDescription()
        {
        }

        virtual std::vector<int> getIndices(size_t aFace) const = 0;
        virtual size_t getIndexCount() const = 0;

        size_t vertexCount;         /**< @brief Number of vertices*/
        size_t faceCount;           /**< @brief Number of faces*/
        const float* vertexArray;   /**< @brief Pointer to the vertex table*/
        const float* normalArray;   /**< @brief Pointer to the normal table (optional)*/

    };

    struct TriangleMeshDescription : public DiscreteGeometryDescription
    {
        TriangleMeshDescription();

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

    /** @brief Represents a container to store visual debugging information collected during a visibility query.*/

    class VisualDebugger
    {
    public:
        VisualDebugger();
        ~VisualDebugger();

        /** @brief Return a mesh containing the removed triangles for visualization purposes*/
        const TriangleMeshDescription& getRemovedTriangles() const;

        /** @brief Return a mesh containing the stabbing lines for visualization purposes*/
        const TriangleMeshDescription& getStabbingLines() const;

        /** @brief Return a mesh containing the extremal stabbing lines for visualization purposes*/
        const TriangleMeshDescription& getExtremalStabbingLines() const;
        
        /** @brief Return a mesh containing the sampling lines for visualization purposes*/
        const TriangleMeshDescription& getSamplingLines() const;

        /** @brief Return a mesh containing the removed edges for visualization purposes*/
        const TriangleMeshDescription& getRemovedEdges() const;


        HelperDebugVisualisation* get()  { return mInternalDebugger; }
        const HelperDebugVisualisation* get() const { return mInternalDebugger; }
    private:
        HelperDebugVisualisation* mInternalDebugger;
    };


    /** @brief Configuration of a visibility query.*/

    struct QueryConfiguration
    {
        /** @brief Arithmetic model presision  */
        enum PrecisionType
        {
            FLOAT,       /**< @brief Single floating point aritmetic*/
            DOUBLE       /**< @brief Double floating point aritmetic*/
#ifdef EXACT_ARITHMETIC  
            , EXACT      /**< @brief Exact arithmetic (based on CGAL and LEDA)*/
#endif
        };

        QueryConfiguration()
        {
            silhouetteOptimization = true;
            hyperSphereNormalization = false;
            performFastTest = false;
            testPolytopeForIntersectionWithOccluder = false;
            representativeLineSampling = true;
            precision = DOUBLE;
            detectApertureOnly = true;
            useEmbree = false;
        };

        QueryConfiguration(const QueryConfiguration& other)
        {
            silhouetteOptimization = other.silhouetteOptimization;
            hyperSphereNormalization = other.hyperSphereNormalization;
            performFastTest = other.performFastTest;
            testPolytopeForIntersectionWithOccluder = other.testPolytopeForIntersectionWithOccluder;
            precision = other.precision;
            representativeLineSampling = other.representativeLineSampling;
            detectApertureOnly = other.detectApertureOnly;
            useEmbree = other.useEmbree;
        };

        bool silhouetteOptimization;                  /**< @brief Use silhouette optimization*/
        bool hyperSphereNormalization;                /**< @brief Normalize the Plucker point (projection on unit hypersphere)*/ 
        bool performFastTest;                         /**< @brief Deprecated */
        bool testPolytopeForIntersectionWithOccluder; /**< @brief Deprecated */
        bool representativeLineSampling;              /**< @brief Use a representative line heuristic to sample the apertures instead of using the extremal stabbing lines*/   
        PrecisionType precision;                      /**< @brief Arithmetic model presicion t*/  
        bool detectApertureOnly;                      /**< @brief Stop the query as soon as a visible line has been found*/  
        bool useEmbree;
    };



    /** @brief Result of a visibility computation */
    enum VisibilityResult
    {
        VISIBLE,  /**< @brief Mutually visible*/
        HIDDEN,   /**< @brief Mutually hidden*/
        UNKNOWN,  /**< @brief Visibility not determined yet*/
        FAILURE   /**< @brief Error occured during computation*/
    };

    /**< @brief Compute if two convex source primitives are mutually visible through the occluders contained in a scene
    @param scene: a scene containing the occluders
    @param vertices0: a pointer to the vertices of the first convex primitive source
    @param numVertices0: the number of vertices of the first convex primitive source
    @param vertices1: a pointer to the vertices of the second convex primitive source
    @param numVertices1: the number of vertices of the second convex primitive source
    @param configuration: configuration parameters of the query  (optional)
    @param debugger: container for debug information during computation (optional)
    @return: the mutual visibility of the two source primitives
    */

    VisibilityResult areVisible(GeometryOccluderSet* scene,
                                const float* vertices0, size_t numVertices0, const float* vertices1, size_t numVertices1,
                                const QueryConfiguration& configuration = QueryConfiguration(), 
                                VisualDebugger* debugger = nullptr);

};

#include "visilib.hpp"

/**@mainpage
# Visilib #

This README would normally document whatever steps are necessary to get your application upand running.

### What is this repository for? ###

The visilib project aims to provide an efficient C++ library for exact visibility computation in 3D.

The main feature of the current version is exact occlusion query between two polygons, a polygonand a segment, or between two segments, and allows to answer the fundamental question "Are these two polygons mutually visible?".The principal application is from - region occlusion cullingand exact Potentially Visible Set computation(PVS).The visibility determination is exactand make no simplification, in opposition to conservative or aggressive visibility methods : visilib performs the complete penumbra fusion of the occluders, by the mean of CSG operations in Plücker space.All the technical details are described in the companion paper "A Low Dimensional Framework for Exact Polygon-to-Polygon Occlusion Queries", presented at EGSR 2005 (cf.download section)

###Features###

- Available exact visibility query :
-polygon to polygon
- segment to segment
- segment to polygon
- point to point
- point to segment
- point to polygon

- robust arithmetic computation via exact artithmetic via CGALand LEDA libraries(multiple precision arithmetic combined with interval arithmetic computations)
- fast ray - triangle intersection based on Intel Embree library
- easy to use : header only library

- Visibility performed via CSG(Constructive Solid Geometry) operations in Plucker space :
-Optimal minimum polytope building algorithm in Plucker space
- Polytope splitting algorithm using only skeleton of polytopes(verticesand edges)
- Silhouette optimization algortihm reducing drastically the number of CSG operations
- Guided aperture samplingand early termination : rays are casted in the visibility apertures left by the previously processed occluders, leading to early termination in case of mutual visibility
- Occluder selection using previous queries
- Computational Geometry predicates in Plucker space

###Applications###
- Potentially Visible Set computation(PVS)
- Sound propagation(edge diffraction)
- Shadow and lighting computation
- Robotic

*/


