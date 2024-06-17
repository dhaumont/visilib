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
#include "demo_helper.h"
#include "math_vector_3.h"
#include "helper_triangle_mesh_container.h"
#include "helper_geometry_scene_reader.h"
#include "helper_synthetic_mesh_builder.h"

using namespace visilib;
using namespace visilibDemo;

void DemoHelper::generatePolygon(std::vector<float>& v, size_t vertexCount, float size, float phi, float scaling)
{
    std::vector<MathVector3f> vertices;

    HelperSyntheticMeshBuilder::generateRegularPolygon(vertices, vertexCount);
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


HelperTriangleMeshContainer* DemoHelper::createScene(int s, float globalScalingFactor)
{
    HelperTriangleMeshContainer* myMeshContainer = new HelperTriangleMeshContainer();

    bool rescale = false;
    if (s == 0)
    {
        if (!DemoHelper::load(myMeshContainer, "..//..//demo//data//sphereWithHoles.obj"))
            return nullptr;
        rescale = true;
    }
    else if (s == 1)
    {
        for (int i = 0; i < 3; i++)
        {
            HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateSlot(0.2f, 0.2f, 0.5f, 0.01f);
            HelperSyntheticMeshBuilder::scale(mesh, 2.0);
            HelperSyntheticMeshBuilder::rotate(mesh, 0, (float)M_PI_2, (float)M_PI);
            HelperSyntheticMeshBuilder::translate(mesh, MathVector3f((float)i / 3.0f, 0.0f, 0.0f));
            myMeshContainer->add(mesh);
            rescale = false;
        }
    }
    else if (s == 2)
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
    }
    else if (s == 3)
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
    }
    else if (s == 4)
    {
        HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateSphere(4);
        HelperSyntheticMeshBuilder::scale(mesh, 0.3f);
        HelperSyntheticMeshBuilder::removeFaces(mesh, 100);

        myMeshContainer->add(mesh);
        rescale = true;
    }
    else if (s == 5)
    {
        HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateSphere(4);
        HelperSyntheticMeshBuilder::addRandomness(mesh, 0.03f);
        HelperSyntheticMeshBuilder::scale(mesh, 0.05f);
        HelperSyntheticMeshBuilder::removeFaces(mesh, 100);

        myMeshContainer->add(mesh);
        rescale = true;
    }
    else if (s == 6)
    {
        HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateCube(1);
        HelperSyntheticMeshBuilder::scale(mesh, 0.3f);
     //   HelperSyntheticMeshBuilder::removeFaces(mesh, 100);

        myMeshContainer->add(mesh);
        rescale = true;
    }
    else if (s == 7)
    {
        HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateCube(4);
        HelperSyntheticMeshBuilder::addRandomness(mesh, 0.03f);
        HelperSyntheticMeshBuilder::removeFaces(mesh, 100);

        HelperSyntheticMeshBuilder::scale(mesh, 0.05f);

        myMeshContainer->add(mesh);
        rescale = true;
    }
    else if (s == 8 || s == 9)
    {
        const size_t count = s == 8 ? 0 : 100;

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

VisibilityExactQueryConfiguration::PrecisionType DemoConfiguration::getPrecisionType() const
{
    if (sampling)
    {
        return VisibilityExactQueryConfiguration::AGGRESSIVE;
    }
#if EXACT_ARITHMETIC     
    if (exactArithmetic)
    {
        return VisibilityExactQueryConfiguration::EXACT;
    }
#endif     
    return VisibilityExactQueryConfiguration::DOUBLE;

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
    std::cout << "[Normalization: " << getStatusString(normalization) << "]" << std::endl;
    std::cout << "[Sampling: " << getStatusString(sampling)<< "]";
    std::cout << "[ExactArithmetic: " << getStatusString(exactArithmetic)<< "]" << std::endl;
    
#if EXACT_ARITHMETIC            
    if (exactArithmetic)
        std::cout << "  [Exact arithmetic: ON]";
    else
#endif
        std::cout << "  [Exact arithmetic: OFF]";
#if EMBREE           
    std::cout << "[Embree:" << getStatusString(embree) << "]" << std::endl;
#endif

}


void DemoConfiguration::writeConfig(const std::string& filename)
{
    std::ofstream output(filename);
    output << "vertexCount0 = " << vertexCount0 << std::endl;
    output << "vertexCount1 = " << vertexCount1 << std::endl;
    output << "silhouetteOptimisation  = " << silhouetteOptimisation << std::endl;
    output << "detectApertureOnly  = " << detectApertureOnly << std::endl;
    output << "normalization  = " << normalization << std::endl;
    output << "scaling  = " << scaling << std::endl;
    output << "phi  = " << phi << std::endl;
    output << "eta = " << eta << std::endl;
    output << "sceneIndex = " << sceneIndex << std::endl;
    output << "globalScaling = " << globalScaling << std::endl;
    output << "exactArithmetic = " << exactArithmetic << std::endl;
    output << "sampling = " << sampling << std::endl;
    output << "tolerance = " << tolerance << std::endl;
    output.close();
}

void DemoConfiguration::readConfig(const std::string& filename)
{
    std::ifstream input(filename);
    if (!input.is_open())
        return;
    while (!input.eof())
    {
        std::vector<std::string> tokens;

        HelperGeometrySceneReader::tokenizeNextLine(input, tokens);
        if (tokens.size() != 3)
            continue;

        if (tokens[0] == "vertexCount0") { vertexCount0 = atoi(tokens[2].c_str()); }
        else if (tokens[0] == "vertexCount1") { vertexCount1 = atoi(tokens[2].c_str()); }
        else if (tokens[0] == "silhouetteOptimisation") { silhouetteOptimisation = atoi(tokens[2].c_str()); }
        else if (tokens[0] == "detectApertureOnly") { detectApertureOnly = atoi(tokens[2].c_str()); }
        else if (tokens[0] == "normalization") { normalization = atoi(tokens[2].c_str()); }
        else if (tokens[0] == "scaling") { scaling = (float)atof(tokens[2].c_str()); }
        else if (tokens[0] == "phi") { phi = (float)atof(tokens[2].c_str()); }
        else if (tokens[0] == "eta") { eta = (float)atof(tokens[2].c_str()); }
        else if (tokens[0] == "sceneIndex") { sceneIndex = atoi(tokens[2].c_str()); }
        else if (tokens[0] == "globalScaling") { globalScaling = (float)atof(tokens[2].c_str()); }
        else if (tokens[0] == "sampling") { sampling = atoi(tokens[2].c_str()); }
        else if (tokens[0] == "tolerance") { tolerance =  (double)atof(tokens[2].c_str()); }
        else { V_ASSERT(0); }
    }
    input.close();
}


void DemoConfiguration::displaySummary()
{
    std::string method = sampling ? "Sampling visibility" : "Exact visibility";
    std::cout << method << " [SceneIndex : " << sceneIndex << ", scaling: " << globalScaling << std::endl;
    std::cout << ", v0: " << vertexCount0 << ", vv1: " << vertexCount1 << "; phi:" << phi;
    std::cout << "; precision: "<< getPrecisionType() <<"; tolerance:" << tolerance << "] ";
}

#if EMBREE
RTCDevice SilhouetteContainerEmbree::mDevice = nullptr;
#endif
