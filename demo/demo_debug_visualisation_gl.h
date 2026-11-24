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

#include "visilib.h"
#include "math_vector_3.h"
#include "helper_triangle_mesh_container.h"

namespace visilibDemo
{
    /* Which of the data a query collected are drawn, each one switched on its own in the interface. */
    struct DemoDisplayOptions
    {
        bool scene = true;                    /* the occluders, wireframe and filled */
        bool removedEdges = false;            /* the occluder edges the solver split with */
        bool removedTriangles = false;        /* the occluder triangles the solver processed */
        bool stabbingLines = false;           /* the stabbing lines of the visible cells */
        bool extremalStabbingLines = false;   /* the extremal stabbing lines of the visible cells */
        bool samplingLines = false;           /* the lines cast to collect the occluders */
        bool shadowDiscontinuities = true;    /* the shadow discontinuities traced on the receiver */
        bool apertures = false;               /* the minimum aperture size, around each source */
        bool cullSourceBackFaces = false;     /* draw the query polygons from their front side only */
    };

    class DemoDebugVisualisationGl
    {
    public:
        static size_t drawAperture(const float size, const visilib::MathVector3f& position );
        static size_t draw(const visilib::GeometryAABB& aabb);
        static size_t drawPolygonEdges(const std::vector<float>& vertices);
        static size_t drawPolygonPoints(const std::vector<float>& vertices);
        static size_t drawPolygon(const std::vector<float>& vertices);
        static size_t drawScene(const visilib::HelperTriangleMeshContainer& aScene);
        static size_t drawMesh(const visilib::GeometryTriangleMeshDescription& aMesh);
        static size_t drawLines(const visilib::GeometryTriangleMeshDescription& info);
        static size_t drawTriangles(const visilib::GeometryTriangleMeshDescription& info);

        static size_t display(visilib::HelperVisualDebugger* debugger, const visilib::HelperTriangleMeshContainer& aScene, const std::vector<float>& v0, const std::vector<float>& v1, visilib::VisibilityResult result, const DemoDisplayOptions& options, float aApertureSize);
    };
}
