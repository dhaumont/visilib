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

#include "visilib.h"
#include "helper_debug_visualisation.h"
#include "visibility_exact_query.h"
#include "geometry_convex_polygon.h"
#include "geometry_occluder_set.h"

using namespace visilib;

inline const TriangleMeshDescription& VisualDebugger::getRemovedTriangles() const
{
    return get()->getRemovedTriangles();
}

inline const TriangleMeshDescription& VisualDebugger::getStabbingLines() const
{
    return get()->getStabbingLines();
}

inline const TriangleMeshDescription& VisualDebugger::getExtremalStabbingLines() const
{
    return get()->getExtremalStabbingLines();
}

inline const TriangleMeshDescription& VisualDebugger::getSamplingLines() const
{
    return get()->getSamplingLines();
}

inline const TriangleMeshDescription& VisualDebugger::getRemovedEdges() const
{
    return get()->getRemovedEdges();
}

inline DiscreteGeometryDescription::DiscreteGeometryDescription()
:  vertexCount(0),
   faceCount(0),
   vertexArray(nullptr),
   normalArray(nullptr)
{
}

inline TriangleMeshDescription::TriangleMeshDescription()
: DiscreteGeometryDescription(),
    indexArray(nullptr)
{
}

inline VisualDebugger::VisualDebugger()
{
    mInternalDebugger = new HelperDebugVisualisation();
}

inline VisualDebugger::~VisualDebugger()
{
    delete mInternalDebugger;
}

inline VisibilityResult visilib::areVisible(GeometryOccluderSet* scene, const float* vertices0, size_t numVertices0, const float* vertices1, size_t numVertices1,
    const QueryConfiguration& configuration, VisualDebugger* debugger)
{
    if (vertices0 == 0 || vertices1 == 0)
    {
        std::cerr << "Error: invalid number of vertices" << std::endl;
        return FAILURE;
    }

    if (scene == nullptr || dynamic_cast<GeometryOccluderSet*>(scene) == nullptr)
    {
        std::cerr << "Error: invalid scene" << std::endl;
        return FAILURE;
    }
    if (vertices0 == nullptr || vertices1 == nullptr)
    {
        std::cerr << "Error: invalid vertex array" << std::endl;
        return FAILURE;
    }

    IVisibilityExactQuery* query = nullptr;

    switch (configuration.precision)
    {
#ifdef EXACT_ARITHMETIC
    case QueryConfiguration::EXACT:
        query = new VisibilityExactQuery<MathPlucker6<exact>, exact>(scene, configuration, MathArithmetic<exact>::Tolerance());
        break;
#endif
        extern double _scaling;
    case QueryConfiguration::DOUBLE:
        query = new VisibilityExactQuery<MathPlucker6<double>, double>(scene, configuration, MathArithmetic<double>::Tolerance());
        break;

    default:

        query = new VisibilityExactQuery<MathPlucker6<float>, float>(scene, configuration, MathArithmetic<float>::Tolerance());
        break;
    }

    query->attachVisualisationDebugger(debugger);

    VisibilityResult result = query->arePolygonsVisible(vertices0, numVertices0, vertices1, numVertices1);

    query->displayStatistic();
    delete query;

    return result;
}
