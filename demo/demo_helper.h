/*
Visilib, an open source library for exact visibility computation.
Copyright(C) 2017 by Denis Haumont

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
#include "helper_triangle_mesh_container.h"

namespace visilibDemo
{
    class DemoHelper
    {
    public:
        static void generatePolygon(std::vector<float>& v, size_t vertexCount, float size, float phi, float scaling);
        static bool load(visilib::HelperTriangleMeshContainer* scene, const std::string& fileName, bool rempveDegeneratedTriangles = false);
    };
}
