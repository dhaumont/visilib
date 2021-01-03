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
#include "geometry_mesh_description.h"
#include "geometry_occluder_set.h"
#include "helper_visual_debugger.h"


namespace visilib
{  
    /** @brief Configuration of a visibility query.*/

    struct VisibilityExactQueryConfiguration
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

        VisibilityExactQueryConfiguration()
        {
            silhouetteOptimization = true;
            hyperSphereNormalization = false;
            representativeLineSampling = true;
            precision = DOUBLE;
            detectApertureOnly = true;
            useEmbree = false;
        };

        VisibilityExactQueryConfiguration(const VisibilityExactQueryConfiguration& other)
        {
            silhouetteOptimization = other.silhouetteOptimization;
            hyperSphereNormalization = other.hyperSphereNormalization;
            precision = other.precision;
            representativeLineSampling = other.representativeLineSampling;
            detectApertureOnly = other.detectApertureOnly;
            useEmbree = other.useEmbree;
        };

        bool silhouetteOptimization;                  /**< @brief Use silhouette optimization*/
        bool hyperSphereNormalization;                /**< @brief Normalize the Plucker point (projection on unit hypersphere)*/ 
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
                                const VisibilityExactQueryConfiguration& configuration = VisibilityExactQueryConfiguration(), 
                                HelperVisualDebugger* debugger = nullptr);

};

#include "visilib.hpp"

/**@mainpage
# Visilib #

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


