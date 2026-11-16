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

#include <cstddef>
#include <unordered_map>

#include "math_vector_3.h"

namespace visilib
{
    /** @brief Provenance of one hyperplane of the Plucker polyhedron: the 3D primitive that carries it.

    The polyhedron stores every constraint line in the same pool of Plucker points, and a facet
    description is just a list of indices into it. The shadow discontinuity tracer needs to know
    which of those lines bound a source polygon, which bound an occluder silhouette, and which
    come from the synthetic edges the builder adds to pad a degenerated source (1 or 2 vertices)
    to a valid polygon: an event vertex built on a synthetic edge is a construction artifact,
    not a shadow boundary. */

    struct HyperplaneEventInfo
    {
        enum Kind
        {
            SourceEdge,       /**< @brief The line supports one edge of one of the two source polygons*/
            SilhouetteEdge    /**< @brief The line supports one edge of an occluder silhouette*/
        };

        HyperplaneEventInfo()
            : mKind(SilhouetteEdge),
            mPolygonIndex(-1),
            mEdgeIndex(0),
            mFaceIndex(0),
            mGeometryId(0),
            mTouchesSynthetic(false)
        {
        }

        Kind mKind;                       /**< @brief Source of the constraint line*/
        int mPolygonIndex;                /**< @brief SourceEdge: which query polygon (0 or 1)*/
        size_t mEdgeIndex;                /**< @brief Index of the edge inside the polygon or the mesh face*/
        size_t mFaceIndex;                /**< @brief SilhouetteEdge: index of the mesh face carrying the edge*/
        size_t mGeometryId;               /**< @brief SilhouetteEdge: identifier of the occluder geometry*/
        bool mTouchesSynthetic;           /**< @brief SourceEdge: the edge uses one of the vertices added by the padding of a degenerated source*/
        MathVector3d mA;                  /**< @brief The two endpoints of the 3D edge carrying the hyperplane*/
        MathVector3d mB;
    };

    /** @brief Side table of the Plucker polyhedron: hyperplane index -> 3D primitive carrying it.

    Populated where the hyperplanes are created (the polytope builder for the source edges, the
    solvers for the silhouette edges), and read by the shadow discontinuity tracer. An index with
    no entry is a line the tracer did not see registered: the tracer still classifies it
    geometrically, it just cannot tell source from occluder, nor real from synthetic. */

    class VisibilityEventRegistry
    {
    public:
        void registerHyperplane(size_t anIndex, const HyperplaneEventInfo& anInfo)
        {
            mEntries[anIndex] = anInfo;
        }

        /** @brief Return the provenance of a hyperplane, or nullptr when the index is unknown*/
        const HyperplaneEventInfo* find(size_t anIndex) const
        {
            auto iter = mEntries.find(anIndex);
            return iter == mEntries.end() ? nullptr : &iter->second;
        }

        void clear()
        {
            mEntries.clear();
        }

        size_t size() const
        {
            return mEntries.size();
        }

    private:
        std::unordered_map<size_t, HyperplaneEventInfo> mEntries;
    };
}
