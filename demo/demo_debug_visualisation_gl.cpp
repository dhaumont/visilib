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

#ifdef _WIN32
#include <windows.h>
#endif

#ifdef USE_GLUT
#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION \
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif
#endif
#include "demo_debug_visualisation_gl.h"
#include <iostream>

using namespace visilib;
using namespace visilibDemo;

size_t DemoDebugVisualisationGl::drawPolygon(const std::vector<float>& vertices)
{
#ifdef USE_GLUT
    glBegin(GL_POLYGON);

    for (size_t i = 0; i < vertices.size(); i += 3)
    {
        glVertex3f(vertices[i], vertices[i + 1], vertices[i + 2]);
    }
    glEnd();
#endif
    return vertices.size() / 3;
}

size_t DemoDebugVisualisationGl::drawPolygonEdges(const std::vector<float> & vertices)
{
    const size_t n = vertices.size();
#ifdef USE_GLUT

    glBegin(GL_POLYGON);

    for (size_t i = 0; i < n; i += 3)
    {
        glVertex3f(vertices[i], vertices[i + 1], vertices[i + 2]);
        glVertex3f(vertices[(i + 3) % n], vertices[(i + 4) % n], vertices[(i + 5) % n]);
    }
    glEnd();
#endif
    return vertices.size() / 3;
}

size_t DemoDebugVisualisationGl::drawPolygonPoints(const std::vector<float> & vertices)
{
    const size_t n = vertices.size();
#ifdef USE_GLUT

    glBegin(GL_POINTS);

    for (size_t i = 0; i < n; i += 3)
    {
        glVertex3f(vertices[i], vertices[i + 1], vertices[i + 2]);
    }
    glEnd();
#endif
    return vertices.size() / 3;
}

size_t DemoDebugVisualisationGl::drawTriangles(const GeometryTriangleMeshDescription & aMesh)
{
#ifdef USE_GLUT

    glBegin(GL_TRIANGLES);

    for (size_t i = 0; i < aMesh.vertexCount * 3; i += 3)
    {
        glVertex3fv(&aMesh.vertexArray[i]);
    }

    glEnd();
#endif
    return aMesh.vertexCount / 3;
}

size_t DemoDebugVisualisationGl::drawLines(const GeometryTriangleMeshDescription & aMesh)
{
#ifdef USE_GLUT

    glBegin(GL_LINES);

    for (size_t i = 0; i < aMesh.vertexCount * 3; i += 3)
    {
        glVertex3fv(&aMesh.vertexArray[i]);
    }

    glEnd();
#endif
    return aMesh.vertexCount / 2;
}

size_t DemoDebugVisualisationGl::drawScene(const HelperTriangleMeshContainer & aScene)
{
    size_t count = 0;

    const std::vector<size_t>& ids = aScene.getGeometryIds();

    for (size_t i = 0; i < ids.size(); i++)
    {
        GeometryDiscreteMeshDescription* mesh = aScene.createTriangleMeshDescription(ids[i]);
        GeometryTriangleMeshDescription* triangleMeshDescription = dynamic_cast<GeometryTriangleMeshDescription*>(mesh);
        if (triangleMeshDescription)
            count += drawMesh(*triangleMeshDescription);
        delete mesh;
    }
    return count;
}

size_t DemoDebugVisualisationGl::drawMesh(const GeometryTriangleMeshDescription & aMesh)
{
#ifdef USE_GLUT

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

    if (aMesh.normalArray != nullptr)
    {
        glNormalPointer(GL_FLOAT, 0, &aMesh.normalArray[0]);
    }

    glVertexPointer(3, GL_FLOAT, 0, &aMesh.vertexArray[0]);

    glDrawElements(GL_TRIANGLES, (int)aMesh.getIndexCount(), GL_UNSIGNED_INT, &aMesh.indexArray[0]);

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
#endif
    return aMesh.faceCount;
}

size_t DemoDebugVisualisationGl::display(HelperVisualDebugger * debugger, const visilib::HelperTriangleMeshContainer & scene, const std::vector<float> & v0, const std::vector<float> & v1, VisibilityResult result, int drawGeometryType)
{
#ifdef USE_GLUT

    glEnable(GL_LIGHTING);

    glEnable(GL_LIGHT0);

    glEnable(GL_DEPTH_TEST);

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glEnable(GL_BLEND);

    glEnable(GL_COLOR_MATERIAL);

    glClearColor(0.1f, 0.0f, 0.0f, 0.0f);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glDisable(GL_LIGHTING);

    glDisable(GL_CULL_FACE);

    glColor3f(0.0f, 1.0f, 1.0f);

    glLineWidth(3);

    DemoDebugVisualisationGl::drawLines(debugger->getRemovedEdges());

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glColor3f(0.9f, 0.9f, 0.9f);

    DemoDebugVisualisationGl::drawTriangles(debugger->getRemovedTriangles());

    if (drawGeometryType < 2)
    {
        glLineWidth(2);

        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        glColor3f(0.6f, 0.6f, 0.6f);

        DemoDebugVisualisationGl::drawScene(scene);

        glDisable(GL_LIGHTING);

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        glColor3f(0.4f, 0.4f, 0.4f);

        DemoDebugVisualisationGl::drawScene(scene);
    }

    switch (result)
    {
    case UNKNOWN:glColor3f(1.0f, 1.0f, 1.0f); break;
    case FAILURE:glColor3f(1.0f, 0.0f, 0.0f); break;
    case VISIBLE:glColor3f(0.0f, 0.9f, 0.7f); break;
    case HIDDEN:glColor3f(0.2f, 0.2f, 0.8f); break;
    }

    if (drawGeometryType < 3)
    {
        glDisable(GL_CULL_FACE);
    }
    else
    { 
       glEnable(GL_CULL_FACE);
    }
    DemoDebugVisualisationGl::drawPolygon(v0);
    DemoDebugVisualisationGl::drawPolygon(v1);

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glLineWidth(3);

    glColor3f(0.9f, 0.9f, 0.9f);

    DemoDebugVisualisationGl::drawPolygonEdges(v0);
    DemoDebugVisualisationGl::drawPolygonEdges(v1);

    glColor3f(255.0f, 255.0f, 0.0f);

    glPointSize(3);
    DemoDebugVisualisationGl::drawPolygonPoints(v0);
    DemoDebugVisualisationGl::drawPolygonPoints(v1);
 
    glLineWidth(1);
    glColor3f(0.0f, 255.0f, 255.0f);
    DemoDebugVisualisationGl::drawLines(debugger->getStabbingLines());
 
    glColor3f(0.0f, 255.0f, 0.0f);
    DemoDebugVisualisationGl::drawLines(debugger->getExtremalStabbingLines());

    glColor3f(0.0f, 0.0f, 255.0f);
    if (drawGeometryType == 1)
    {
        DemoDebugVisualisationGl::drawLines(debugger->getSamplingLines());
    }
#endif
    return 0;
}