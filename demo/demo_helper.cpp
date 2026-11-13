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
        "NOISY_CUBE", "CUBES_AND_GRID", "CUBES_AND_GRID_100"
    };

    const size_t count = sizeof(names) / sizeof(names[0]);
    return scene >= 0 && static_cast<size_t>(scene) < count ? names[scene] : "UNKNOWN_SCENE";
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
