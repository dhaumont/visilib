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
#define _USE_MATH_DEFINES
#include <math.h>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <algorithm>

#ifdef _WIN32
#include <windows.h>
#include <direct.h>
#else
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif
#include "demo_helper.h"
#include "math_vector_3.h"
#include "math_arithmetic.h"
#include "helper_triangle_mesh_container.h"
#include "helper_geometry_scene_reader.h"
#include "helper_synthetic_mesh_builder.h"
#include "visilib_core.h"
/* Included by path relative to this file so that the test target, which also compiles
demo_helper.cpp without adding demo/external/json to its include dirs, still finds it. */
#include "external/json/json.hpp"

using namespace visilib;
using namespace visilibDemo;

void DemoHelper::generatePolygon(std::vector<float>& v, size_t vertexCount, float size, float phi, float scaling)
{
    std::vector<MathVector3f> vertices;
    std::vector<int> indices;

    HelperSyntheticMeshBuilder::generateRegularPolygon(indices, vertices, vertexCount);
    HelperSyntheticMeshBuilder::scale(vertices, size * scaling);
    HelperSyntheticMeshBuilder::rotate(vertices, 0, (float)M_PI_2, phi);
    HelperSyntheticMeshBuilder::translate(vertices,  MathVector3f(scaling * cos(phi), scaling * sin(phi), 0));

    v.clear();
    for (size_t i = 0; i < vertices.size(); i++)
    {
        v.push_back(vertices[i].x); v.push_back(vertices[i].y); v.push_back(vertices[i].z);
    }
}

bool DemoHelper::load(visilib::HelperTriangleMeshContainer* scene, const std::string& fileName, bool removeDegeneratedTriangles)
{
    HelperGeometrySceneReader reader(scene);

    bool result = reader.readFileObj(fileName);

    if (removeDegeneratedTriangles)
    {
        scene->removeDegeneratedTriangles(0.000001f);
    }
    return result;
}

int DemoHelper::getSubdivision(DemoConfiguration::SCENE_TYPE s)
{
    const double size = 0.07;
    switch(s)
    {
        case DemoConfiguration::SINGLE_TRIANGLE_01: return 0;
        case DemoConfiguration::SINGLE_TRIANGLE_02: return 1;
        case DemoConfiguration::SINGLE_TRIANGLE_03: return 2;
        case DemoConfiguration::SINGLE_TRIANGLE_04: return 3;        
    }

    return 1.0;
}

double DemoHelper::getApertureSize(DemoConfiguration::SCENE_TYPE s)
{
    const double size = 0.07;
    switch(s)
    {
        case DemoConfiguration::SLOT_OFF_AXIS_01: return size;
        case DemoConfiguration::SLOT_OFF_AXIS_02: return size/2;
        case DemoConfiguration::SLOT_OFF_AXIS_03: return size/4;
        case DemoConfiguration::SLOT_OFF_AXIS_04: return size/8;
        case DemoConfiguration::SLOT_OFF_AXIS_05: return size/16;
    }

    return 1.0;
}

namespace
{

    void addFusionLayers(HelperTriangleMeshContainer* container, float width, float height, float slit,
                         float stagger, float slide)
    {
        const int layerCount = 10;
        const float pitch = height + slit;
        const float firstDepth = -0.45f;
        const float depthStep = 0.1f;

        for (int i = 0; i < layerCount; i++)
        {
            /* The alternating half pitch puts the plates of one layer opposite the slits of the
            layer turned a quarter of a turn from it; the stagger does the same for the layers of
            the one orientation, which are four turns apart. */
            const float base = 0.5f * pitch * (i % 2 == 0 ? -1.0f : 1.0f)
                + stagger * ((i / 2) % 2);
            const float turn = (float)M_PI_2 * (float)i;
            const float along = slide * ((float)i - 0.5f * (float)(layerCount - 1));

            for (int k = 0; k < 5; k++)
            {
                /* The bar runs along y with its thickness along x, and the turn the other facing
                plates use puts it in a plane x = constant: width across, height up. The plate is set
                in its row along that height, the whole row is slipped along its own length and swung
                about the view axis, and only then is it moved out to the depth of its layer. */
                HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateBar(width, height);
                HelperSyntheticMeshBuilder::rotate(mesh, 0, (float)M_PI_2, (float)M_PI);
                HelperSyntheticMeshBuilder::translate(mesh, MathVector3f(0.0f, along, base + ((float)k - 2.0f) * pitch));
                HelperSyntheticMeshBuilder::rotate(mesh, turn, 0.0f, 0.0f);
                HelperSyntheticMeshBuilder::translate(mesh, MathVector3f(firstDepth + (float)i * depthStep, 0.0f, 0.0f));

                container->add(mesh);
            }
        }
    }
}

HelperTriangleMeshContainer* DemoHelper::createScene(DemoConfiguration::SCENE_TYPE s, float globalScalingFactor)
{
    HelperTriangleMeshContainer* myMeshContainer = new HelperTriangleMeshContainer();

    bool rescale = false;

    V_TRACE("INIT:" << s);
    switch(s)
    {
     //   case DemoConfiguration::SPHERE_WITH_HOLES:
       // {
         //   if (!DemoHelper::load(myMeshContainer, "..//..//demo//data//sphereWithHoles.obj"))
           //     return nullptr;
           // rescale = true;
        //}break;

        case DemoConfiguration::SLOT_OFF_AXIS_01:
        case DemoConfiguration::SLOT_OFF_AXIS_02:
        case DemoConfiguration::SLOT_OFF_AXIS_03:
        case DemoConfiguration::SLOT_OFF_AXIS_04:
        case DemoConfiguration::SLOT_OFF_AXIS_05:
        {
            for (int i = 0; i < 2; i++)
            {
                double apertureSize = getApertureSize(s);
                HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateSlot(0.0f, 0.0f, apertureSize, apertureSize);
                HelperSyntheticMeshBuilder::scale(mesh, 2.0);
                HelperSyntheticMeshBuilder::rotate(mesh, 0, (float)M_PI_2, (float)M_PI);
                HelperSyntheticMeshBuilder::translate(mesh, MathVector3f((float)i / 3.0f, 0.0f, 0.0f));
                myMeshContainer->add(mesh);
                rescale = false;
            }
        }break;
       case  DemoConfiguration::SLOT_ON_AXIS_01:
        {
            const size_t count = 2;
            for (int i = 0; i < count; i++)
            {
                HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateSlot(0.0f, 0.0f, 0.03f, 0.03f);
                HelperSyntheticMeshBuilder::scale(mesh, 2.0);
                HelperSyntheticMeshBuilder::rotate(mesh, 0.0, (float)M_PI_2, (float)M_PI);

                if (count > 1)
                    HelperSyntheticMeshBuilder::translate(mesh, MathVector3f(-0.5f + (float)i / (float)(count + 1), 0, 0));
                myMeshContainer->add(mesh);
                rescale = false;
            }
        }break;
       case  DemoConfiguration::REGULAR_GRID:
        {
            const size_t count = 1;
            for (int i = 0; i < count; i++)
            {
                //   GeometryMesh* mesh = HelperSyntheticMeshBuilder::generateSlot(0.0,0.0,0.8,0.02);
                HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateRegularGrid(0);

                HelperSyntheticMeshBuilder::rotate(mesh, 0.0, (float)M_PI_2, (float)M_PI);
                HelperSyntheticMeshBuilder::rotate(mesh, (float)i * (float)M_PI / (5.f * (float)count), 0.0, 0.0);
                HelperSyntheticMeshBuilder::scale(mesh, 2.0);

                HelperSyntheticMeshBuilder::translate(mesh, MathVector3f(0.0f, 0.0f, -1.0f));
                if (count > 1)
                    HelperSyntheticMeshBuilder::translate(mesh, MathVector3f(-0.5f + (float)i / (float)(count + 1), 0.0f, 0.0f));
                myMeshContainer->add(mesh);
                rescale = false;
            }
        }break;
       case DemoConfiguration::SPHERE_WITH_HOLES:
        {
            HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateSphere(4);
            HelperSyntheticMeshBuilder::scale(mesh, 0.3f);
            HelperSyntheticMeshBuilder::removeFaces(mesh, 100);

            myMeshContainer->add(mesh);
            rescale = true;
        }break;
       case DemoConfiguration::NOISY_SPHERE:
        {
            HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateSphere(4);
            HelperSyntheticMeshBuilder::addRandomness(mesh, 0.03f);
            HelperSyntheticMeshBuilder::scale(mesh, 0.05f);
            HelperSyntheticMeshBuilder::removeFaces(mesh, 100);

            myMeshContainer->add(mesh);
            rescale = true;
        }break;
       case  DemoConfiguration::SINGLE_TRIANGLE_01:
       case  DemoConfiguration::SINGLE_TRIANGLE_02:
       case  DemoConfiguration::SINGLE_TRIANGLE_03:
       case  DemoConfiguration::SINGLE_TRIANGLE_04:
        {
            int subdivision = getSubdivision(s);
            HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateRegularPolygon(3., subdivision);
            HelperSyntheticMeshBuilder::rotate(mesh, 0, (float)M_PI_2, (float)M_PI);
     
            HelperSyntheticMeshBuilder::scale(mesh, 0.1);
        
            myMeshContainer->add(mesh);
            rescale = false;
        }break;

       case  DemoConfiguration::SIMPLE_CUBE:
        {
            HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateCube(1);
            HelperSyntheticMeshBuilder::scale(mesh, 0.3f);
        //   HelperSyntheticMeshBuilder::removeFaces(mesh, 100);

            myMeshContainer->add(mesh);
            rescale = true;
        }break;
       case  DemoConfiguration::NOISY_CUBE:
        {
            HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateCube(4);
            HelperSyntheticMeshBuilder::addRandomness(mesh, 0.03f);
            HelperSyntheticMeshBuilder::removeFaces(mesh, 100);

            HelperSyntheticMeshBuilder::scale(mesh, 0.05f);

            myMeshContainer->add(mesh);
            rescale = true;
        }break;
       case  DemoConfiguration::CUBES_AND_GRID:
       case DemoConfiguration::CUBES_AND_GRID_100:
        {
            const size_t count = s == DemoConfiguration::CUBES_AND_GRID ? 0 : 100;

            for (int i = 0; i < count; i++)
            {
                //   GeometryMesh* mesh = HelperSyntheticMeshBuilder::generateSlot(0.0,0.0,0.8,0.02);
                HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateCube(2);
                HelperSyntheticMeshBuilder::scale(mesh, 0.1f);

                HelperSyntheticMeshBuilder::addRandomness(mesh, 0.012f);
                //         HelperSyntheticMeshBuilder::removeFaces(mesh, 100);

                HelperSyntheticMeshBuilder::rotate(mesh, 0.0, (float)M_PI_2, (float)M_PI);
                MathVector3f random(-0.5f + (float)rand() / (float)(RAND_MAX), -0.5f + (float)rand() / (float)(RAND_MAX), -0.5f + (float)rand() / (float)(RAND_MAX));
                //   random *= 0.3;
                HelperSyntheticMeshBuilder::translate(mesh, random);

                myMeshContainer->add(mesh);

            }
            HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateRegularGrid(6);

            HelperSyntheticMeshBuilder::rotate(mesh, 0.0, (float)M_PI_2, (float)M_PI);
            //         HelperSyntheticMeshBuilder::addRandomness(mesh, 0.02);
        //   HelperSyntheticMeshBuilder::removeFaces(mesh, 100);
            HelperSyntheticMeshBuilder::removeFaces(mesh, 100);

            myMeshContainer->add(mesh);


            rescale = false;
        }break;

        case DemoConfiguration::DEGENERATE_FAN_PLATE:
        {
            /* A facing plate cut as a fan, so every spoke ends on the apex, and the ring edges
            come in pairs parallel to each other and lying in the one plane of the plate. The apex
            is set at the centre of the view axis: any sight line through it satisfies all the
            spoke constraints at once, which is the plane of the lines through a point, a plane
            lying on the Klein quadric itself. The facets of the polytope are then not in general
            position, and the cut applied near the apex has to be read against a facet that is
            already there. */
            HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateRegularPolygon(6, 2);
            HelperSyntheticMeshBuilder::scale(mesh, 0.3f);
            HelperSyntheticMeshBuilder::rotate(mesh, 0, (float)M_PI_2, (float)M_PI);

            const MathVector3f apex = mesh->getVertices()[0];
            HelperSyntheticMeshBuilder::translate(mesh, MathVector3f(-apex.x, -apex.y, -apex.z));

            myMeshContainer->add(mesh);
            rescale = false;
        }break;

        case DemoConfiguration::DEGENERATE_HAIRLINE_SLOT:
        {
            /* Two facing plates, each pierced by a crack two thousandths across, the second slid
            sideways by half a crack so that the two cracks only partly line up. The two edges
            bounding one crack are parallel lines a hair apart in the same plane, so their
            hyperplanes are as close to duplicates as the arrangement allows without being equal,
            and the aperture they leave is far below the minimum aperture size: the tolerance alone
            decides what the answer is. */
            const float crack = 0.002f;

            for (int i = 0; i < 2; i++)
            {
                HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateSlot(0.0f, 0.0f, crack, crack);
                HelperSyntheticMeshBuilder::scale(mesh, 2.0f);
                HelperSyntheticMeshBuilder::rotate(mesh, 0, (float)M_PI_2, (float)M_PI);
                HelperSyntheticMeshBuilder::translate(mesh, MathVector3f(0.33f * (float)i, 0.5f * crack * (float)i, 0.0f));
                myMeshContainer->add(mesh);
            }

            rescale = false;
        }break;

        case DemoConfiguration::DEGENERATE_ALIGNED_STACK:
        {
            /* Two plates hinged on the view axis itself, one a hair behind the other, opened like
            a partly shut book across the way the two query polygons face each other. The edges
            they present to the sight are the same line to within two thousandths, so the second
            plate cuts the polytope by a hyperplane repeating a facet it already carries, and the
            sight line down the middle of the axis is on the boundary of the answer rather than
            inside it. */
            const float offset = 0.002f;
            const float hinge[2] = { 0.35f, (float)M_PI - 0.35f };

            for (int i = 0; i < 2; i++)
            {
                HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateRegularGrid(0);
                HelperSyntheticMeshBuilder::scale(mesh, 0.5f);
                /* One edge of the plate is laid along the local x axis, so that the rotation
                about x hinges the plate on that edge and the edge stays on the view axis. */
                HelperSyntheticMeshBuilder::translate(mesh, MathVector3f(0.0f, 0.25f, 0.0f));
                HelperSyntheticMeshBuilder::rotate(mesh, hinge[i], 0.0f, 0.0f);
                HelperSyntheticMeshBuilder::translate(mesh, MathVector3f(0.15f * (float)i, 0.0f, offset * (float)i));
                myMeshContainer->add(mesh);
            }

            rescale = false;
        }break;

        case DemoConfiguration::DEGENERATE_REGULUS_FRAME:
        {
            /* Six slats, each of them a plate parallel to the two query polygons, set at six
            depths across the way between them and turned in their own plane by an angle that grows
            with the depth, so that the six stripes they cast run in six slightly different
            directions across the sight region. Their long axes are six generators of one ruling of
            the saddle z = x * y, and the lines of one ruling have the whole other ruling common to
            them, every one of those transversals running along the view axis, so the six
            constraints leave the search a flat lying on the quadric instead of pinning down a
            vertex. Lifting each slat two thousandths off the saddle, alternately up and down, takes
            the flat away and leaves the splitter a pair of extremal stabbing lines to tell apart.
            Being plates rather than blades seen edge on, they are wide enough to cover the sight
            region and block most of what goes through. */
            const float slat[6] = { -0.5f, -0.3f, -0.1f, 0.1f, 0.3f, 0.5f };
            const float lift = 0.002f;

            for (int i = 0; i < 6; i++)
            {
                HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateBar(1.6f, 0.16f);
                /* The same turn the other facing plates use, which puts the bar in a plane
                x = constant with its long side horizontal, and then a spin about that plane's
                normal, atan(t) for the slat at depth t, tilting its long axis to (0, 1, t). */
                HelperSyntheticMeshBuilder::rotate(mesh, 0, (float)M_PI_2, (float)M_PI);
                HelperSyntheticMeshBuilder::rotate(mesh, atanf(slat[i]), 0.0f, 0.0f);
                HelperSyntheticMeshBuilder::translate(mesh, MathVector3f(slat[i], 0.0f, (i % 2 == 0 ? lift : -lift)));
                myMeshContainer->add(mesh);
            }

            rescale = false;
        }break;

        case DemoConfiguration::DEGENERATE_CUBE_ON_GRID:
        {
            /* A grid plate, whose interior edges are coplanar with the plate and parallel two by
            two, and a cube standing on it with one edge a thousandth aside one of those interior
            lines. Two occluders therefore contribute facets whose hyperplanes nearly coincide,
            which is the seam a real scene has wherever two meshes meet, and the silhouette filter
            keeps both because they belong to different meshes. */
            HelperTriangleMesh* plate = HelperSyntheticMeshBuilder::generateRegularGrid(1);
            HelperSyntheticMeshBuilder::scale(plate, 0.6f);
            HelperSyntheticMeshBuilder::rotate(plate, 0, (float)M_PI_2, (float)M_PI);
            myMeshContainer->add(plate);

            HelperTriangleMesh* cube = HelperSyntheticMeshBuilder::generateCube(0);
            HelperSyntheticMeshBuilder::scale(cube, 0.3f);
            /* The edges along z then lie at x and y of -0.15 and 0.15, so half a side brings two
            of them into the plane of the plate; the odd thousandth keeps one of them just off the
            plate's own interior edge instead of exactly on it, which is the near miss the search
            cannot resolve. */
            HelperSyntheticMeshBuilder::translate(cube, MathVector3f(0.15f, 0.151f, 0.0f));
            myMeshContainer->add(cube);

            rescale = false;
        }break;

        case DemoConfiguration::DEGENERATE_TRIANGLE_FAN_DISK:
        {
            /* A disk cut as a fan whose shared vertex is its center, one layer of elongated
            triangles, each a slender wedge from the center out to the rim. Twenty four of them all
            meet at the one point, and that point sits on the view axis where the sight lines are
            densest: any sight line through the center answers all the spoke constraints at once,
            which is the plane of the lines through a point, a plane lying on the Klein quadric
            itself. The rim edges face each other across the center, parallel and in the one plane
            of the disk. */
            HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateTriangleFan(24, 0);
            HelperSyntheticMeshBuilder::scale(mesh, 0.35f);
            HelperSyntheticMeshBuilder::rotate(mesh, 0, (float)M_PI_2, (float)M_PI);

            myMeshContainer->add(mesh);
            rescale = false;
        }break;

        case DemoConfiguration::DEGENERATE_CUBE_GRID:
        {
            /* Five by five by five cubes on a regular lattice, every face axis aligned, two
            hundredths apart, so that the only way through is the lattice of gaps in three
            dimensions. Each face is parallel to twenty of its neighbours and its edges are collinear
            with theirs, so the arrangement carries whole families of facets whose hyperplanes agree
            to within the width of a gap, and a sight line threaded through one gap passes faces that
            cannot be told apart. */
            const int count = 5;
            const float side = 0.1f;
            const float pitch = 0.12f;

            for (int i = 0; i < count; i++)
            {
                for (int j = 0; j < count; j++)
                {
                    for (int k = 0; k < count; k++)
                    {
                        HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateCube(0);
                        HelperSyntheticMeshBuilder::scale(mesh, side);
                        HelperSyntheticMeshBuilder::translate(mesh,
                            MathVector3f(((float)i - (float)(count - 1) * 0.5f) * pitch,
                                         ((float)j - (float)(count - 1) * 0.5f) * pitch,
                                         ((float)k - (float)(count - 1) * 0.5f) * pitch));

                        myMeshContainer->add(mesh);
                    }
                }
            }

            rescale = false;
        }break;

        case DemoConfiguration::DEGENERATE_CRACKED_CUBES:
        {
            /* Two cubes of one size, all faces parallel to each other, four thousandths apart so
            that the two faces they turn to each other open a crack through the middle of the sight.
            Unlike a hole in a plate this is a slit with depth: a sight line has to hold its two
            thousandths of y across the whole width of the cubes to get through, and the two faces
            are the same plane to within the crack, so the second cube cuts the polytope with a
            hyperplane repeating a facet the first one already supplied. */
            const float crack = 0.004f;
            const float side = 0.4f;

            for (int i = 0; i < 2; i++)
            {
                HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateCube(0);
                HelperSyntheticMeshBuilder::scale(mesh, side);
                HelperSyntheticMeshBuilder::translate(mesh,
                    MathVector3f(0.0f, (i == 0 ? -1.0f : 1.0f) * (side * 0.5f + crack * 0.5f), 0.0f));

                myMeshContainer->add(mesh);
            }

            rescale = false;
        }break;

        case DemoConfiguration::DEGENERATE_OCCLUSION_FUSION:
        {
            /* Ten layers of five small plates each, every plate parallel to the two query polygons:
            each of the fifty plates blocks some lines, no single layer blocks the query, and it is
            the fusion of them that comes close to leaving nothing through, which is the answer no
            occluder can be asked for on its own. Being parallel plates stacked in rows, their edges
            run in families of parallel lines whose hyperplanes stay close to one another over the
            whole search. */
            addFusionLayers(myMeshContainer, 0.5f, 0.048f, 0.04f, 0.0f, 0.02f);

            rescale = false;
        }break;

        case DemoConfiguration::DEGENERATE_OCCLUSION_FUSION_BLOCKED:
        {
            /* The same ten turned layers, with the plates widened from half to nine tenths across
            and from a fiftyth to a seventh and a half high, the slits between them narrowed from
            four hundredths to two, and every second layer of an orientation given half a pitch of
            stagger. The rows then reach past the sight region in every direction they are turned
            to, and a ray that holds a slit at one depth has to jump six hundredths of height to
            hold one two tenths of depth further on, which is more than the two query polygons a
            fifth apart and four tenths of height wide allow it. The fusion closes, and the query is
            hidden by the pile rather than by any one plate in it. */
            addFusionLayers(myMeshContainer, 0.9f, 0.14f, 0.02f, 0.08f, 0.025f);

            rescale = false;
        }break;

        case DemoConfiguration::DEGENERATE_SPIKE:
        {
            /* A ball with thirty thorns stuck in it, a small cactus, centered on the view axis. Every
            thorn is a thin tetrahedron whose apex is the center of the ball, so the three long
            edges of all of them run through that one point and the point is where the sight lines
            are densest: the facets they cut are the facets of the lines through a point, a plane
            lying on the Klein quadric, and the sight line down the axis itself, which goes through
            that point, answers all of them at once. The thorns stand out in five rings of six, two
            rings about each pole and one about the equator, so that the polar ones point almost down
            the axis and their tips nearly touch it, and each ring is turned half a step from the one
            before. A thorn is four thousandths thick at its base for a tenth of length, and shorter
            than the radius of the ball it grows from, so it reads as a thorn rather than as an arm:
            its three edges are a cluster of nearly identical facets, and the ninety of them a family
            that no tolerance sets apart from another. */
            const float ballRadius = 0.15f;
            const float needleLength = 0.10f;
            const float needleRadius = 0.004f;
            const float rings[5] = { 0.3f, 0.75f, (float)M_PI_2, (float)M_PI - 0.75f, (float)M_PI - 0.3f };

            HelperTriangleMesh* ball = HelperSyntheticMeshBuilder::generateSphere(1);
            HelperSyntheticMeshBuilder::scale(ball, ballRadius);
            myMeshContainer->add(ball);

            for (int ring = 0; ring < 5; ring++)
            {
                for (int j = 0; j < 6; j++)
                {
                    const float azimuth = 2.f * (float)M_PI * (float)j / 6.f
                        + (ring % 2 == 0 ? 0.0f : (float)M_PI / 6.f);

                    /* Turning the needle about x by the polar angle tips its axis out of z, and the
                    roll afterwards swings that axis around the view axis, which is how the ring is
                    completed. */
                    HelperTriangleMesh* needle = HelperSyntheticMeshBuilder::generateSpike(needleLength, needleRadius);
                    HelperSyntheticMeshBuilder::rotate(needle, rings[ring], 0.0f, azimuth);

                    myMeshContainer->add(needle);
                }
            }

            rescale = false;
        }break;

        case DemoConfiguration::DEGENERATE_SPIRAL:
        {
            /* The regulus frame carried all the way round: twelve plates parallel to the two query
            polygons, at twelve depths across the corridor, each one turned in its own plane a
            thirtieth of a full turn further than the last, so that the stripe it casts has swept
            thirty degrees on from the previous and the twelve together have gone the whole way
            round. Their long axes are twelve generators of one ruling of the saddle z = x * y,
            lifted two thousandths off it one plate up and two down, so the flat of common
            transversals the ruling leaves is only just taken away, and the stripes cover the sight
            region in every direction at once. */
            const float plate[12] = { -0.55f, -0.45f, -0.35f, -0.25f, -0.15f, -0.05f,
                                       0.05f,  0.15f,  0.25f,  0.35f,  0.45f,  0.55f };
            const float lift = 0.002f;

            for (int i = 0; i < 12; i++)
            {
                HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateBar(1.6f, 0.16f);
                HelperSyntheticMeshBuilder::rotate(mesh, 0, (float)M_PI_2, (float)M_PI);
                HelperSyntheticMeshBuilder::rotate(mesh, 2.f * (float)M_PI * (float)i / 12.f, 0.0f, 0.0f);
                HelperSyntheticMeshBuilder::translate(mesh, MathVector3f(plate[i], 0.0f, (i % 2 == 0 ? lift : -lift)));

                myMeshContainer->add(mesh);
            }

            rescale = false;
        }break;

        case DemoConfiguration::MODEL_ABACUS:
        case DemoConfiguration::MODEL_BIRCH01:
        case DemoConfiguration::MODEL_CAVATAP2:
        case DemoConfiguration::MODEL_CENTERPLANT:
        case DemoConfiguration::MODEL_CENTERPLANT3:
        case DemoConfiguration::MODEL_CERCLE:
        case DemoConfiguration::MODEL_CHAINLNK:
        case DemoConfiguration::MODEL_CONE:
        case DemoConfiguration::MODEL_DNA:
        case DemoConfiguration::MODEL_GLASSES:
        case DemoConfiguration::MODEL_HAND:
        case DemoConfiguration::MODEL_SMALL:
        case DemoConfiguration::MODEL_SPHERE_WITH_HOLES:
        case DemoConfiguration::MODEL_SPHERE_WITH_HOLES2:
        case DemoConfiguration::MODEL_TENNIS_RACKET:
        {
            /* The models shipped in demo/data, each one read as it stands and brought to the unit
            box by the rescale below, so that whatever units a file was written in, the thing set
            between the two query polygons is one box wide. Faces left degenerate by the reading are
            dropped on the way in, which is what the hand built scenes never need and real data
            often does. A file that is not there gives no scene at all rather than an empty box. */
            if (!DemoHelper::load(myMeshContainer, DemoHelper::dataFolder() + "/" + DemoConfiguration::modelName(s), true))
                return nullptr;

            rescale = true;
        }break;
        }
    if (rescale)
    {
        HelperSyntheticMeshBuilder::rescaleToUnitBox(myMeshContainer);
    }
    HelperSyntheticMeshBuilder::scale(myMeshContainer, globalScalingFactor);

    return myMeshContainer;
}

GeometryOccluderSet* DemoHelper::createOccluderSet(HelperTriangleMeshContainer* aContainer)
{
    GeometryOccluderSet* occluderSet = new GeometryOccluderSet();
    for (size_t index = 0; index < aContainer->getGeometryCount(); index++)
    {
        GeometryDiscreteMeshDescription* info = aContainer->createTriangleMeshDescription(index);
        occluderSet->addOccluder(info);
    }
    occluderSet->prepare();

    return occluderSet;
}

void DemoHelper::exportQueryToObj(const std::string& fileName, const std::vector<float>& v0, const std::vector<float>& v1, const HelperTriangleMeshContainer & aScene, float aScaling)
{
    std::ofstream myOutput(fileName.c_str());
    auto myMeshArray  = aScene.getMeshArray();
    int myOffset = 1;

    DemoHelper::appendPolygonRepresentationToFileObj(myOutput, myOffset, v0, 0.01 * aScaling);
    DemoHelper::appendPolygonRepresentationToFileObj(myOutput, myOffset, v1, 0.01 * aScaling);

    for (auto iter = myMeshArray.begin(); iter != myMeshArray.end(); iter++)
    {
        HelperGeometrySceneReader::appendMeshToFileObj(myOutput, myOffset, (*iter)->getVertices(), (*iter)->getIndices());
    }

    myOutput.close();   
}

void DemoHelper::appendPolygonRepresentationToFileObj(std::ofstream& stream, int& anOffset, const std::vector<float>& vertices, float aScaling)
{
    size_t myVertexCount = vertices.size() / 3;
    if (myVertexCount == 0)
        return;
    size_t i = 1;
    if (myVertexCount == 1)
    {
            //export a sphere to represent single points
        MathVector3f translation(vertices[0], vertices[1], vertices[2]);
        std::vector<int> sphere_indices;
        std::vector<MathVector3f> sphere_vertices;
        HelperSyntheticMeshBuilder::generateSphere(2, sphere_indices,  sphere_vertices);
        HelperSyntheticMeshBuilder::scale(sphere_vertices, aScaling);
        HelperSyntheticMeshBuilder::translate(sphere_vertices, translation);
        
        HelperGeometrySceneReader::appendMeshToFileObj(stream, anOffset, sphere_vertices, sphere_indices);            
    }
    else if (myVertexCount == 2)
    {
        std::vector<int> myIndices = {0,1};
        std::vector<MathVector3f> myVertices = 
        {
            MathVector3f(vertices[0], vertices[1], vertices[2]),
            MathVector3f(vertices[3], vertices[4], vertices[5])
        };

        HelperGeometrySceneReader::appendSegmentsToFileObj(stream, anOffset, myVertices, myIndices);      
    }
    else 
    {
        HelperGeometrySceneReader::appendPolygonToFileObj(stream, anOffset, vertices);
    }
}


void DemoHelper::configureDemoConfiguration(const std::string& name, DemoConfiguration& configuration)
{
}

void DemoConfiguration::displaySettings()
{
    std::cout << std::endl << "Current Demo Settings: " << std::endl;

    std::cout << "  [OccluderSet index:" << sceneIndex << "]";
    std::cout << "[Source Size: " << scaling << "]";
    std::cout << "[Source Vertices: " << vertexCount0 << " and " << vertexCount1 << "]";
    std::cout << "[GlobalScaling:" << globalScaling << "]";
    std::cout << "[Tolerance:" << tolerance << "]";
    std::cout << std::endl;

    std::cout << "  [Early stop: " << getStatusString(detectApertureOnly) << "]";
    std::cout << "[Silhouette: " << getStatusString(silhouetteOptimisation) << "]";
    std::cout << "[Normalization: " << getStatusString(normalization) << "]";
    std::cout << "[Rep. line: " << (optimisedRepresentativeLine ? "OPTIMISED" : "GRAVITY CENTER") << "]";
    std::cout << "[SolverType: " << toStr(solverType)<< "]";
  
    std::cout << "  [Arithmetic: " << toStr(precisionType) << "]" << std::endl;
#if EMBREE           
    std::cout << "[Embree:" << getStatusString(embree) << "]" << std::endl;
#endif

}


namespace
{
    /* Keys written by writeConfig: used to tell a typo apart from a key that belongs to
    a newer build of the demo. */
    const char* const cConfigKeys[] = {
        "vertexCount0", "vertexCount1", "silhouetteOptimisation", "detectApertureOnly",
        "normalization", "scaling", "phi", "eta", "sceneIndex", "globalScaling",
        "precisionType", "solverType", "tolerance", "minimumApertureSize", "confidenceValue",
        "expected"
#if EMBREE
        , "embree"
#endif
    };

    bool isKnownConfigKey(const std::string& key)
    {
        for (size_t i = 0; i < sizeof(cConfigKeys) / sizeof(cConfigKeys[0]); ++i)
        {
            if (key == cConfigKeys[i])
                return true;
        }
        return false;
    }

    /* Widening a float to double prints artifacts (0.21f becomes 0.20999999344348907).
    Emit the shortest decimal that reads back as the very same float, so the file stays
    readable without losing precision. */
    double asJsonNumber(float value)
    {
        for (int digits = 1; digits <= std::numeric_limits<float>::max_digits10; ++digits)
        {
            std::ostringstream text;
            text << std::setprecision(digits) << value;
            if (std::strtof(text.str().c_str(), nullptr) != value)
                continue;
            return std::strtod(text.str().c_str(), nullptr);
        }
        return static_cast<double>(value);
    }

    /* A present key carrying the wrong JSON type is reported and ignored, so that a
    hand edited file cannot silently half-apply. */
    void warnUnexpectedType(const nlohmann::json& root, const char* key)
    {
        if (root.contains(key))
            std::cerr << "config: key '" << key << "' has an unexpected type (" << root[key].type_name() << ")" << std::endl;
    }

    void readBool(const nlohmann::json& root, const char* key, bool& value)
    {
        if (root.contains(key) && root[key].is_boolean()) value = root[key].get<bool>();
        else warnUnexpectedType(root, key);
    }

    void readString(const nlohmann::json& root, const char* key, std::string& value)
    {
        if (root.contains(key) && root[key].is_string()) value = root[key].get<std::string>();
        else warnUnexpectedType(root, key);
    }

    /* Numbers are read through double so that a hand edited file cannot overflow the
    target type: non finite values, or values outside the target range, are reported and
    dropped instead of quietly turning into inf or a wrapped integer. */
    template <typename T>
    void readNumber(const nlohmann::json& root, const char* key, T& value)
    {
        if (!root.contains(key))
            return;
        if (root[key].is_number())
        {
            const double candidate = root[key].get<double>();
            if (std::isfinite(candidate)
                && candidate <= static_cast<double>(std::numeric_limits<T>::max())
                && candidate >= static_cast<double>(std::numeric_limits<T>::lowest()))
            {
                value = static_cast<T>(candidate);
                return;
            }
            std::cerr << "config: key '" << key << "' is not a finite value in range, keeping "
                << value << std::endl;
            return;
        }
        warnUnexpectedType(root, key);
    }

    /* Vertex counts must be a whole number in a range the demo can actually build:
    accepting -5 as size_t would wrap around to 1.8e19 and hang the mesh generator. */
    const size_t cMaxVertexCount = 1024;

    void readCount(const nlohmann::json& root, const char* key, size_t& value)
    {
        if (!root.contains(key))
            return;
        if (root[key].is_number_integer())
        {
            const long long count = root[key].get<long long>();
            if (count >= 1 && count <= static_cast<long long>(cMaxVertexCount))
            {
                value = static_cast<size_t>(count);
                return;
            }
        }
        std::cerr << "config: key '" << key << "' must be an integer in [1, "
            << cMaxVertexCount << "], keeping " << value << std::endl;
    }

    /* Enum fields are stored as their numeric value, as in the previous text format, and
    must stay inside the enumerators the scene builders know about. */
    template <typename T>
    void readEnum(const nlohmann::json& root, const char* key, T& value, int count)
    {
        if (!root.contains(key))
            return;
        if (root[key].is_number_integer())
        {
            const long long index = root[key].get<long long>();
            if (index >= 0 && index < static_cast<long long>(count))
            {
                value = static_cast<T>(index);
                return;
            }
        }
        std::cerr << "config: key '" << key << "' must be an integer in [0, " << count
            << "[, keeping " << static_cast<int>(value) << std::endl;
    }

    std::string tolowerCopy(const std::string& text)
    {
        std::string result = text;
        for (size_t i = 0; i < result.size(); ++i)
        {
            if (result[i] >= 'A' && result[i] <= 'Z') result[i] = char(result[i] - 'A' + 'a');
        }
        return result;
    }

    /* File names are built from numbers, and a dot in a file name reads like an
    extension: 0.21 is written 0_21. Lower case, because two names differing only by
    case would be a mess on a case insensitive filesystem. */
    std::string fileNameTag(double value)
    {
        std::ostringstream text;
        text << value;

        std::string tag = tolowerCopy(text.str());
        for (size_t i = 0; i < tag.size(); ++i)
        {
            if (tag[i] == '.' || tag[i] == '-' || tag[i] == '+') tag[i] = '_';
        }
        return tag;
    }

    /* toStr spells the enums out for the console ("double 64 bits; tol = 1e-11"), which
    makes unreadable and space ridden file names, so the names use short codes. */
    std::string solverTag(VisibilityExactQueryConfiguration::SolverType solver)
    {
        const std::string name = tolowerCopy(DemoConfiguration::toStr(solver));

        if (name.find("aperture") != std::string::npos) return "finder";
        if (name.find("sequential") != std::string::npos) return "sequential";
        if (name.find("monte") != std::string::npos) return "mc";
        return "solver";
    }

    bool readFile(const std::string& path, std::string& content)
    {
        std::ifstream input(path.c_str());
        if (!input.is_open())
            return false;

        std::ostringstream buffer;
        buffer << input.rdbuf();
        content = buffer.str();
        return true;
    }

    /* Only the last component may be missing, and an existing folder is not an error. */
    void ensureFolder(const std::string& folder)
    {
#ifdef _WIN32
        _mkdir(folder.c_str());
#else
        mkdir(folder.c_str(), 0777);
#endif
    }

    /* Trailing blanks carry no meaning in a JSON document, but they do differ when two
    files are compared as text. */
    std::string trimmed(const std::string& text)
    {
        size_t end = text.size();
        while (end > 0 && (text[end - 1] == '\n' || text[end - 1] == '\r' || text[end - 1] == ' ' || text[end - 1] == '\t'))
            --end;
        return text.substr(0, end);
    }
}

std::string DemoConfiguration::jsonString() const
{
    nlohmann::json root;

    root["vertexCount0"] = vertexCount0;
    root["vertexCount1"] = vertexCount1;
    root["silhouetteOptimisation"] = silhouetteOptimisation;
    root["detectApertureOnly"] = detectApertureOnly;
    root["normalization"] = normalization;
    root["scaling"] = asJsonNumber(scaling);
    root["phi"] = asJsonNumber(phi);
    root["eta"] = asJsonNumber(eta);
    root["sceneIndex"] = sceneIndex;
    root["globalScaling"] = asJsonNumber(globalScaling);
    root["precisionType"] = precisionType;
    root["solverType"] = solverType;
    root["tolerance"] = tolerance;
    root["minimumApertureSize"] = minimumApertureSize;
    root["confidenceValue"] = confidenceValue;

    /* Only written when set: a file exported straight from the interface carries no
    expectation, and its absence must not show up as an empty key. */
    if (!expected.empty())
        root["expected"] = expected;
#if EMBREE
    root["embree"] = embree;
#endif

    return root.dump(2);
}

void DemoConfiguration::writeConfig(const std::string& filename)
{
    std::ofstream output(filename);
    if (!output.is_open())
    {
        std::cerr << "config: cannot write '" << filename << "'" << std::endl;
        return;
    }
    output << jsonString() << std::endl;
    output.close();
}

const char* DemoConfiguration::sceneName(SCENE_TYPE scene)
{
    /* Same order as the SCENE_TYPE enumerators and as the Scene combo of the interface. */
    static const char* names[] =
    {
        "SINGLE_TRIANGLE_01", "SINGLE_TRIANGLE_02", "SINGLE_TRIANGLE_03", "SINGLE_TRIANGLE_04",
        "SLOT_OFF_AXIS_01", "SLOT_OFF_AXIS_02", "SLOT_OFF_AXIS_03", "SLOT_OFF_AXIS_04", "SLOT_OFF_AXIS_05",
        "SLOT_ON_AXIS_01", "SPHERE_WITH_HOLES", "REGULAR_GRID", "NOISY_SPHERE", "SIMPLE_CUBE",
        "NOISY_CUBE", "CUBES_AND_GRID", "CUBES_AND_GRID_100",
        "DEGENERATE_FAN_PLATE", "DEGENERATE_HAIRLINE_SLOT", "DEGENERATE_ALIGNED_STACK",
        "DEGENERATE_REGULUS_FRAME", "DEGENERATE_CUBE_ON_GRID",
        "DEGENERATE_TRIANGLE_FAN_DISK", "DEGENERATE_CUBE_GRID", "DEGENERATE_CRACKED_CUBES",
        "DEGENERATE_OCCLUSION_FUSION", "DEGENERATE_SPIKE", "DEGENERATE_SPIRAL",
        "DEGENERATE_OCCLUSION_FUSION_BLOCKED",
        "MODEL_ABACUS", "MODEL_BIRCH01", "MODEL_CAVATAP2", "MODEL_CENTERPLANT", "MODEL_CENTERPLANT3",
        "MODEL_CERCLE", "MODEL_CHAINLNK", "MODEL_CONE", "MODEL_DNA", "MODEL_GLASSES", "MODEL_HAND",
        "MODEL_SMALL", "MODEL_SPHERE_WITH_HOLES", "MODEL_SPHERE_WITH_HOLES2", "MODEL_TENNIS_RACKET"
    };

    const size_t count = sizeof(names) / sizeof(names[0]);
    static_assert(count == DemoConfiguration::COUNT, "one name per scene, in the order of SCENE_TYPE");
    return scene >= 0 && static_cast<size_t>(scene) < count ? names[scene] : "UNKNOWN_SCENE";
}

const char* DemoConfiguration::modelName(SCENE_TYPE scene)
{
    /* The files of demo/data, from the first MODEL_* scene on, in the order of that block. The
    names are the ones the folder holds, capitalisation and all, so that a scene can be matched to
    its file at a glance. */
    static const char* files[] =
    {
        "abacus.obj", "BIRCH01.obj", "CAVATAP2.obj", "centerplant.obj", "centerplant3.obj",
        "cercle.obj", "CHAINLNK.obj", "cone.obj", "DNA.obj", "GLASSES.obj", "hand.obj",
        "small.obj", "sphereWithHoles.obj", "sphereWithHoles2.obj", "tennis-racket.obj"
    };

    const size_t first = MODEL_ABACUS;
    const size_t count = sizeof(files) / sizeof(files[0]);
    static_assert(static_cast<size_t>(MODEL_TENNIS_RACKET - MODEL_ABACUS) + 1 == count, "one file per model scene");

    return scene >= 0 && static_cast<size_t>(scene) >= first
        && static_cast<size_t>(scene) < first + count ? files[static_cast<size_t>(scene) - first] : "";
}

std::string DemoConfiguration::caseFileName() const
{
    std::ostringstream file;

    file << tolowerCopy(sceneName(sceneIndex))
        << "__s" << vertexCount0 << "_" << vertexCount1
        << "__sc" << fileNameTag(scaling)
        << "__phi" << fileNameTag(phi)
        << "__eta" << fileNameTag(eta)
        << "__gs" << fileNameTag(globalScaling)
        << "__f" << (silhouetteOptimisation ? "1" : "0")
                 << (normalization ? "1" : "0")
                 << (detectApertureOnly ? "1" : "0")
        << "__" << DemoHelper::precisionTag(precisionType)
        << "__" << solverTag(solverType)
        << "__tol" << fileNameTag(tolerance)
        << "__ap" << fileNameTag(minimumApertureSize)
        << "__cf" << fileNameTag(confidenceValue);
#if EMBREE
    file << "__emb" << (embree ? "1" : "0");
#endif
    file << ".json";

    return file.str();
}

bool DemoConfiguration::writeTestCase(const std::string& folder, std::string& message) const
{
    const std::string name = caseFileName();
    const std::string path = folder + "/" + name;
    const std::string body = jsonString();

    std::string existing;
    if (readFile(path, existing))
    {
        message = trimmed(existing) == body
            ? name + " already holds these settings"
            : name + " already exists with different content, keeping it as it is";
        return false;
    }

    ensureFolder(folder);

    std::ofstream output(path.c_str());
    if (!output.is_open())
    {
        message = "could not write " + path;
        return false;
    }

    output << body << std::endl;
    output.close();

    message = "wrote " + path;
    return true;
}

void DemoConfiguration::readConfig(const std::string& filename)
{
    std::ifstream input(filename);
    if (!input.is_open())
        return;

    /* Never throw on a hand edited file: is_discarded reports the syntax error and the
    current configuration is kept untouched. */
    nlohmann::json root = nlohmann::json::parse(input, nullptr, false);
    input.close();

    if (root.is_discarded() || !root.is_object())
    {
        std::cerr << "config: '" << filename << "' is not a valid JSON object" << std::endl;
        return;
    }

    readCount(root, "vertexCount0", vertexCount0);
    readCount(root, "vertexCount1", vertexCount1);
    readBool(root, "silhouetteOptimisation", silhouetteOptimisation);
    readBool(root, "detectApertureOnly", detectApertureOnly);
    readBool(root, "normalization", normalization);
    readNumber(root, "scaling", scaling);
    readNumber(root, "phi", phi);
    readNumber(root, "eta", eta);
    readEnum(root, "sceneIndex", sceneIndex, DemoConfiguration::COUNT);
    readNumber(root, "globalScaling", globalScaling);
    readEnum(root, "precisionType", precisionType, VisibilityExactQueryConfiguration::PRECISION_COUNT);
    readEnum(root, "solverType", solverType, VisibilityExactQueryConfiguration::SOLVER_COUNT);
    readNumber(root, "tolerance", tolerance);
    readNumber(root, "minimumApertureSize", minimumApertureSize);
    readNumber(root, "confidenceValue", confidenceValue);
    readString(root, "expected", expected);
#if EMBREE
    readBool(root, "embree", embree);
#endif

    /* Unknown keys are reported and skipped rather than rejected: a file written by a
    newer build must still load in an older one, and vice versa. */
    for (auto it = root.begin(); it != root.end(); ++it)
    {
        if (!isKnownConfigKey(it.key()))
            std::cout << "config: ignoring unknown key '" << it.key() << "'" << std::endl;
    }
}

const std::string DemoConfiguration::toStr(VisibilityExactQueryConfiguration::PrecisionType p)
{
    std::stringstream stream;
    switch (p) 
    {
    case VisibilityExactQueryConfiguration::FLOAT:
        stream <<  "FLOAT "  << MathArithmetic<float>::bitsCount()<< " bits; tol = " << MathArithmetic<float>::Tolerance();
        break;
    case VisibilityExactQueryConfiguration::DOUBLE:
         stream << "DOUBLE " << MathArithmetic<double>::bitsCount()<< " bits; tol = " << MathArithmetic<double>::Tolerance();
        break;
#ifdef ENABLE_LEDA
    case VisibilityExactQueryConfiguration::LEDA_REAL:
        stream << "LEDA_REAL " <<MathArithmetic<MathLedaReal>::bitsCount()<< " bits; tol = " << MathArithmetic<MathLedaReal>::Tolerance();
        break;
#endif
#ifdef ENABLE_GMP
    case VisibilityExactQueryConfiguration::GMP_FLOAT:
        stream << "GMP_FLOAT "  <<MathArithmetic<MathGmpFloat>::bitsCount()<< " bits; tol = " << MathArithmetic<MathGmpFloat>::Tolerance();
        break;
    case VisibilityExactQueryConfiguration::GMP_RATIONAL:
        stream << "GMP_RATIONAL; tol = " << MathArithmetic<MathGmpRational>::Tolerance();
        break;

#endif
#ifdef ENABLE_MPFR
    case VisibilityExactQueryConfiguration::MPFR:
        stream << "MPFR " << MathArithmetic<MathMpfr>::bitsCount() << " bits; tol = " << MathArithmetic<MathMpfr>::Tolerance();
        break;

#endif
    default:
        stream << "UNKNOWN";
        break;
    }

    return stream.str();
}


const std::string DemoConfiguration::toStr(VisibilityExactQueryConfiguration::SolverType p)
{
    switch (p) 
    {
    case VisibilityExactQueryConfiguration::EXACT_APERTURE_FINDER:
        return "EXACT_APERTURE_FINDER";
        break;
    case VisibilityExactQueryConfiguration::EXACT_SEQUENTIAL_SOLVER:
        return "EXACT_SEQUENTIAL_SOLVER";
        break;
    case VisibilityExactQueryConfiguration::MONTE_CARLO:
        return "MONTE_CARLO";
        break;
    }
    return "UNKNOWN";
}

const std::string DemoConfiguration::toStr(visilib::VisibilityResult result)
{
    switch (result)
    {
        case visilib::VISIBLE: return "VISIBLE";
        case visilib::HIDDEN:  return "HIDDEN";
        case visilib::UNKNOWN: return "UNKNOWN";
        default:               return "FAILURE";
    }
}

/* Anything that is not VISIBLE or HIDDEN reads back as UNKNOWN, which the
consistency test takes for "no pin": a hand edited typo stays honest instead of
failing with a wrong expectation. */
visilib::VisibilityResult DemoConfiguration::resultFromName(const std::string& name)
{
    if (name == "VISIBLE") return visilib::VISIBLE;
    if (name == "HIDDEN")  return visilib::HIDDEN;
    return visilib::UNKNOWN;
}

void DemoConfiguration::displaySummary()
{
    std::string method = toStr(solverType);
    std::cout << method << " [SceneIndex : " << sceneIndex << ", scaling: " << globalScaling << std::endl;
    std::cout << ", v0: " << vertexCount0 << ", vv1: " << vertexCount1 << "; phi:" << phi;
    std::cout << "; precision: "<< toStr(precisionType) <<"; tolerance:" << tolerance << "] ";
}

#if EMBREE
RTCDevice SilhouetteContainerEmbree::mDevice = nullptr;
#endif

/* toStr spells the enums out for the console ("double 64 bits; tol = 1e-11"), which
makes unreadable and space ridden file names, so the names use short codes. */
std::string DemoHelper::precisionTag(VisibilityExactQueryConfiguration::PrecisionType precision)
{
    const std::string name = tolowerCopy(DemoConfiguration::toStr(precision));

    if (name.find("leda") != std::string::npos) return "LEDA";
    if (name.find("rational") != std::string::npos) return "GMPQ";
    if (name.find("gmp") != std::string::npos) return "GMPF";
    if (name.find("mpfr") != std::string::npos) return "MPFR";
    if (name.find("float") != std::string::npos) return "SP";
    return "DP";
}

std::string DemoHelper::testCaseFolder()
{
    /* The demo runs from wherever it is started, so the folder it saves into has to be
    known at configure time rather than guessed from the current directory. */
#ifdef TEST_CASES_DIR
    return TEST_CASES_DIR;
#else
    return "test_cases";
#endif
}

std::string DemoHelper::dataFolder()
{
#ifdef DEMO_DATA_DIR
    return DEMO_DATA_DIR;
#else
    return "demo/data";
#endif
}

std::vector<std::string> DemoHelper::listTestCaseFiles(const std::string& folder)
{
    std::vector<std::string> files;

#ifdef _WIN32
    WIN32_FIND_DATAA entry;
    HANDLE search = FindFirstFileA((folder + "\*.json").c_str(), &entry);
    if (search == INVALID_HANDLE_VALUE)
        return files;

    do
    {
        files.push_back(entry.cFileName);
    } while (FindNextFileA(search, &entry) != 0);
    FindClose(search);
#else
    DIR* directory = opendir(folder.c_str());
    if (directory == nullptr)
        return files;

    /* readdir is only asked about names, the .json suffix is the filter. */
    struct dirent* entry = nullptr;
    while ((entry = readdir(directory)) != nullptr)
    {
        const std::string name(entry->d_name);
        if (name.size() > 5 && name.compare(name.size() - 5, 5, ".json") == 0)
            files.push_back(name);
    }
    closedir(directory);
#endif

    /* The list is what the interface shows, so it must come out in the same order every
    time; readdir has no ordering promise. */
    std::sort(files.begin(), files.end());
    return files;
}
