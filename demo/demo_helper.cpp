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

#include "demo_helper.h"
#include "math_vector_3.h"
#include "helper_triangle_mesh_container.h"
#include "helper_geometry_scene_reader.h"
#include "helper_synthetic_mesh_builder.h"

using namespace visilib;
using namespace visilibDemo;

void DemoHelper::generatePolygon(std::vector<float>& v, size_t vertexCount, float size, float phi, float scaling)
{
    const float PI = 3.14159265358979323846f;
    std::vector<MathVector3f> vertices;

    HelperSyntheticMeshBuilder::generateRegularPolygon(vertices, vertexCount);
    HelperSyntheticMeshBuilder::scale(vertices, size * scaling);
    HelperSyntheticMeshBuilder::rotate(vertices, 0, PI / 2, phi);
    HelperSyntheticMeshBuilder::translate(vertices,  MathVector3f(scaling * cos(phi), scaling * sin(phi), 0));

    v.clear();
    for (size_t i = 0; i < vertices.size(); i++)
    {
        v.push_back(vertices[i].x); v.push_back(vertices[i].y); v.push_back(vertices[i].z);
    }
}

bool DemoHelper::load(HelperTriangleMeshContainer* scene, const std::string& fileName, bool removeDegeneratedTriangles)
{
    HelperGeometrySceneReader reader(scene);

    bool result = reader.readFileObj(fileName);

    if (removeDegeneratedTriangles)
    {
        scene->removeDegeneratedTriangles(0.000001f);
    }
    return result;
}
