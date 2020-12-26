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

#include <vector>

#include "math_vector_3.h"
#include "visilib.h"
#include <fstream>
#include "helper_debug_visualisation.h"

namespace visilib
{

    /** @brief Represents a container to store visual debugging information collected during a visibility query.*/

    class HelperVisualDebugger
    {
    public:
        HelperVisualDebugger();
        ~HelperVisualDebugger();

        /** @brief Return a mesh containing the removed triangles for visualization purposes*/
        const GeometryTriangleMeshDescription& getRemovedTriangles() const;

        /** @brief Return a mesh containing the stabbing lines for visualization purposes*/
        const GeometryTriangleMeshDescription& getStabbingLines() const;

        /** @brief Return a mesh containing the extremal stabbing lines for visualization purposes*/
        const GeometryTriangleMeshDescription& getExtremalStabbingLines() const;

        /** @brief Return a mesh containing the sampling lines for visualization purposes*/
        const GeometryTriangleMeshDescription& getSamplingLines() const;

        /** @brief Return a mesh containing the removed edges for visualization purposes*/
        const GeometryTriangleMeshDescription& getRemovedEdges() const;


        HelperDebugVisualisation* get() { return mInternalDebugger; }
        const HelperDebugVisualisation* get() const { return mInternalDebugger; }
    private:
        HelperDebugVisualisation* mInternalDebugger;
    };


    inline const GeometryTriangleMeshDescription& HelperVisualDebugger::getRemovedTriangles() const
    {
        return get()->getRemovedTriangles();
    }

    inline const GeometryTriangleMeshDescription& HelperVisualDebugger::getStabbingLines() const
    {
        return get()->getStabbingLines();
    }

    inline const GeometryTriangleMeshDescription& HelperVisualDebugger::getExtremalStabbingLines() const
    {
        return get()->getExtremalStabbingLines();
    }

    inline const GeometryTriangleMeshDescription& HelperVisualDebugger::getSamplingLines() const
    {
        return get()->getSamplingLines();
    }

    inline const GeometryTriangleMeshDescription& HelperVisualDebugger::getRemovedEdges() const
    {
        return get()->getRemovedEdges();
    }


    inline HelperVisualDebugger::HelperVisualDebugger()
    {
        mInternalDebugger = new HelperDebugVisualisation();
    }

    inline HelperVisualDebugger::~HelperVisualDebugger()
    {
        delete mInternalDebugger;
    }

}