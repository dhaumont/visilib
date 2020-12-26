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

#ifdef _WIN32
#include <windows.h>
#endif

#include <GL/gl.h>
#include <GL/glut.h>

#include <string>
#include <iostream>
#include <fstream>

#include "xmmintrin.h"
#include "pmmintrin.h"

#include "helper_triangle_mesh_container.h"
#include "silhouette_container_embree.h"
#include "helper_synthetic_mesh_builder.h"
#include "helper_geometry_scene_reader.h"
#include "demo_viewer_glut.h"
#include "demo_debug_visualisation_gl.h"
#include "demo_helper.h"

using namespace visilib;
using namespace std;
using namespace visilibDemo;

namespace visilibDemo
{
    class VisilibDemoMain;

    class VisilibDemoMain
    {
    public:
        VisilibDemoMain()
        {
        }

        std::string getStatusString(bool enable)
        {
            return enable ? "ON" : "OFF";
        }

        bool init()
        {

            _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
            _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
            readConfig("config.txt");
            forceDisplay = true;

            gluLookAt(1, 1, 1,
                0, 0, 0,
                1, 0, 0);

            return initScene(sceneIndex);
        }
        bool initScene(int s)
        {
            delete debugger;
            delete occluderSet;
            delete meshContainer;
            srand(0);
            debugger = new VisualDebugger();

            meshContainer = new HelperTriangleMeshContainer();

            bool rescale = false;
            if (s == 0)
            {
                if (!DemoHelper::load(meshContainer, "..//..//demo//data//sphereWithHoles.obj"))
                    return false;
                rescale = true;
            }
            else if (s == 1)
            {
                for (int i = 0; i < 3; i++)
                {
                    HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateSlot(0.2f, 0.2f, 0.5f, 0.01f);
                    HelperSyntheticMeshBuilder::scale(mesh, 2.0);
                    HelperSyntheticMeshBuilder::rotate(mesh, 0, 3.14f / 2, 3.14f);
                    HelperSyntheticMeshBuilder::translate(mesh, MathVector3f((float)i / 3.0f, 0.0f, 0.0f));
                    meshContainer->add(mesh);
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
                    HelperSyntheticMeshBuilder::rotate(mesh, 0.0, 3.14f / 2, 3.14f);

                    if (count > 1)
                        HelperSyntheticMeshBuilder::translate(mesh, MathVector3f(-0.5f + (float)i / (float)(count + 1), 0, 0));
                    meshContainer->add(mesh);
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

                    HelperSyntheticMeshBuilder::rotate(mesh, 0.0, 3.14f / 2, 3.14f);
                    HelperSyntheticMeshBuilder::rotate(mesh, (float)i * 3.14f / (5 * (float)count), 0.0, 0.0);
                    HelperSyntheticMeshBuilder::scale(mesh, 2.0);

                    HelperSyntheticMeshBuilder::translate(mesh, MathVector3f(0.0f, 0.0f, -1.0f));
                    if (count > 1)
                        HelperSyntheticMeshBuilder::translate(mesh, MathVector3f(-0.5f + (float)i / (float)(count + 1), 0.0f, 0.0f));
                    meshContainer->add(mesh);
                    rescale = false;
                }
            }
            else if (s == 4)
            {
                HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateSphere(4);
                HelperSyntheticMeshBuilder::scale(mesh, 0.3f);
                HelperSyntheticMeshBuilder::removeFaces(mesh, 100);

                meshContainer->add(mesh);
                rescale = true;
            }
            else if (s == 5)
            {
                HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateSphere(4);
                HelperSyntheticMeshBuilder::addRandomness(mesh, 0.03f);
                HelperSyntheticMeshBuilder::scale(mesh, 0.05f);
                HelperSyntheticMeshBuilder::removeFaces(mesh, 100);

                meshContainer->add(mesh);
                rescale = true;
            }
            else if (s == 6)
            {
                HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateCube(4);
                HelperSyntheticMeshBuilder::scale(mesh, 0.3f);
                HelperSyntheticMeshBuilder::removeFaces(mesh, 100);

                meshContainer->add(mesh);
                rescale = true;
            }
            else if (s == 7)
            {
                HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateCube(4);
                HelperSyntheticMeshBuilder::addRandomness(mesh, 0.03f);
                HelperSyntheticMeshBuilder::removeFaces(mesh, 100);

                HelperSyntheticMeshBuilder::scale(mesh, 0.05f);

                meshContainer->add(mesh);
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

                    HelperSyntheticMeshBuilder::rotate(mesh, 0.0, 3.14f / 2, 3.14f);
                    MathVector3f random(-0.5f + (float)rand() / (float)(RAND_MAX), -0.5f + (float)rand() / (float)(RAND_MAX), -0.5f + (float)rand() / (float)(RAND_MAX));
                    //   random *= 0.3;
                    HelperSyntheticMeshBuilder::translate(mesh, random);

                    meshContainer->add(mesh);

                }
                HelperTriangleMesh* mesh = HelperSyntheticMeshBuilder::generateRegularGrid(6);

                HelperSyntheticMeshBuilder::rotate(mesh, 0.0, 3.14f / 2, 3.14f);
                //         HelperSyntheticMeshBuilder::addRandomness(mesh, 0.02);
             //   HelperSyntheticMeshBuilder::removeFaces(mesh, 100);
                HelperSyntheticMeshBuilder::removeFaces(mesh, 100);

                meshContainer->add(mesh);


                rescale = false;
            }

            if (rescale)
            {
                HelperSyntheticMeshBuilder::rescaleToUnitBox(meshContainer);
            }
            HelperSyntheticMeshBuilder::scale(meshContainer, globalScaling);

            occluderSet = new GeometryOccluderSet();
            for (size_t index = 0; index < meshContainer->getGeometryCount(); index++)
            {
                DiscreteGeometryDescription* info = meshContainer->createTriangleMeshDescription(index);
                occluderSet->addOccluder(info);
            }
            occluderSet->prepare();
            return true;
        }

        void resolveVisibility()
        {
            QueryConfiguration config;
            config.silhouetteOptimization = silhouetteOptimisation;
            config.hyperSphereNormalization = normalization;
            config.precision = precisionType;
            config.representativeLineSampling = representativeLineSampling;
            config.detectApertureOnly = detectApertureOnly;
#if EMBREE 
            config.useEmbree = embree;
#endif
            result = visilib::areVisible(occluderSet, &v0[0], v0.size() / 3, &v1[0], v1.size() / 3, config, debugger);
        }

        void animate()
        {
            if (animated)
            {
                phi += 0.005f;
                eta += 0.001f;
                forceDisplay = true;
            }

            if (forceDisplay)
            {
                DemoHelper::generatePolygon(v0, vertexCount, scaling, phi - 3.14519f, globalScaling);
                DemoHelper::generatePolygon(v1, 3, scaling, phi, globalScaling);

                resolveVisibility();

                forceDisplay = false;
            }

            glutPostRedisplay();
        }

        void writeConfig(const std::string & filename)
        {
            std::ofstream output(filename);
            output << "vertexCount = " << vertexCount << std::endl;
            output << "silhouetteOptimisation  = " << silhouetteOptimisation << std::endl;
            output << "detectApertureOnly  = " << detectApertureOnly << std::endl;
            output << "drawGeometryType  = " << drawGeometryType << std::endl;
            output << "representativeLineSampling  = " << representativeLineSampling << std::endl;
            output << "animated  = " << animated << std::endl;
            output << "normalization  = " << normalization << std::endl;
            output << "scaling  = " << scaling << std::endl;
            output << "phi  = " << phi << std::endl;
            output << "eta = " << eta << std::endl;
            output << "sceneIndex = " << sceneIndex << std::endl;
            output << "globalScaling = " << globalScaling << std::endl;
            output << "precisionType = " << precisionType << std::endl;

            output.close();
        }

        void readConfig(const std::string & filename)
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

                if (tokens[0] == "vertexCount") { vertexCount = atoi(tokens[2].c_str()); }
                else if (tokens[0] == "silhouetteOptimisation") { silhouetteOptimisation = atoi(tokens[2].c_str()); }
                else if (tokens[0] == "detectApertureOnly") { detectApertureOnly = atoi(tokens[2].c_str()); }
                else if (tokens[0] == "drawGeometryType") { drawGeometryType = atoi(tokens[2].c_str()); }
                else if (tokens[0] == "representativeLineSampling") { representativeLineSampling = atoi(tokens[2].c_str()); }
                else if (tokens[0] == "animated") { animated = atoi(tokens[2].c_str()); }
                else if (tokens[0] == "normalization") { normalization = atoi(tokens[2].c_str()); }
                else if (tokens[0] == "scaling") { scaling = (float)atof(tokens[2].c_str()); }
                else if (tokens[0] == "phi") { phi = (float)atof(tokens[2].c_str()); }
                else if (tokens[0] == "eta") { eta = (float)atof(tokens[2].c_str()); }
                else if (tokens[0] == "sceneIndex") { sceneIndex = atoi(tokens[2].c_str()); }
                else if (tokens[0] == "globalScaling") { globalScaling = (float)atof(tokens[2].c_str()); }
                else if (tokens[0] == "precisionType") { precisionType = (QueryConfiguration::PrecisionType)atoi(tokens[2].c_str()); }

                else { V_ASSERT(0); }
            }
            input.close();
        }

        void display()
        {
            DemoDebugVisualisationGl::display(debugger, *meshContainer, v0, v1, result, drawGeometryType);
        }

        void displaySettings()
        {
            std::cout << std::endl << "Current Demo Settings: " << std::endl;

            std::cout << "  [OccluderSet index:" << sceneIndex << "]";
            std::cout << "[Source Size: " << scaling << "]";
            std::cout << "[Source Vertices: " << vertexCount << "]";
            std::cout << "[GlobalScaling:" << globalScaling << "]";
            std::cout << std::endl;

            std::cout << "  [Early stop: " << getStatusString(detectApertureOnly) << "]";
            std::cout << "[Silhouette: " << getStatusString(silhouetteOptimisation) << "]";
            std::cout << "[Middle line: " << getStatusString(representativeLineSampling) << "]";
            std::cout << "[Normalization: " << getStatusString(normalization) << "]" << std::endl;

#if EXACT_ARITHMETIC            
            if (precisionType == QueryConfiguration::EXACT)
                std::cout << "  [Exact arithmetic: ON]";
            else
#endif
                std::cout << "  [Exact arithmetic: OFF]";
#if EMBREE           
            std::cout << "[Embree:" << getStatusString(embree) << "]" << std::endl;
#endif
            
        }

        void writeHelp()
        {
            std::cout << "Visilib 1.0. Demo application " << std::endl;

            std::cout << "  s: enable/disable silhouette optimisation" << std::endl;
            std::cout << "  n: enable/disable nomalization" << std::endl;
#if EXACT_ARITHMETIC
            std::cout << "  e: enable/disable exact arithmetic" << std::endl;
#endif
#if EMBREE
            std::cout << "  g: enable/disable embree ray tracing" << std::endl;
#endif
            std::cout << "  r: enable/disable representative line sampling strategy" << std::endl;
            std::cout << "  f: enable/disable detect aperture only" << std::endl;

            std::cout << "  n: enable/disable Plucker normalization" << std::endl;
            std::cout << "  f: enable/disable fast silhouette rejection test" << std::endl;

            std::cout << "  x: change scene " << std::endl;
            std::cout << "  +/-: increase/decrease scaling of query polygons" << std::endl;
            std::cout << "  1/2: increase/decrease number of vertices of query polygons" << std::endl;

            std::cout << "  w: write config" << std::endl;
            std::cout << "  o: open config" << std::endl;

            std::cout << "  space: start/pause animation" << std::endl;
            std::cout << "  Enter: show/hide geometry" << std::endl;

            std::cout << "  h: write this help" << std::endl;
        }

        void keyboard(unsigned char key, int, int)
        {
            switch (key)
            {
            case 27:  // The escape key
            case 'Q':
            case 'q':
                exit(0);   // Simply exit
                break;

            case '2':
                if (vertexCount < 12)
                    vertexCount++;

                forceDisplay = true;
                break;
            case '1':
                if (vertexCount > 1)
                    vertexCount--;

                forceDisplay = true;
                break;

            case '+':
                if (scaling < 1.00f)
                    scaling += 0.01f;
                forceDisplay = true;
                break;

            case '-':
                if (scaling > 0.02f)
                    scaling -= 0.01f;

                forceDisplay = true;
                break;

            case '*':
                globalScaling *= 2;
                forceDisplay = true;
                setViewPortScaling(globalScaling);
                initScene(sceneIndex);

                break;

            case '/':
                globalScaling /= 2;

                forceDisplay = true;
                setViewPortScaling(globalScaling);
                initScene(sceneIndex);

                break;

            case 'h':

                writeHelp();

                displaySettings();
                break;
            case 's':
                silhouetteOptimisation = !silhouetteOptimisation;
                forceDisplay = true;
                break;
            case 'f':
                detectApertureOnly = !detectApertureOnly;
                forceDisplay = true;

                break;

            case 'x':
                sceneIndex++;
                if (sceneIndex > 9)
                    sceneIndex = 0;
                initScene(sceneIndex);
                forceDisplay = true;
                break;
            case 'r':
                representativeLineSampling = !representativeLineSampling;
                forceDisplay = true;
                break;


#ifdef EXACT_ARITHMETIC
            case 'e':
                precisionType = precisionType == QueryConfiguration::DOUBLE ? QueryConfiguration::EXACT : QueryConfiguration::DOUBLE;
                forceDisplay = true;
#endif
                break;
#if EMBREE
            case 'g':
                embree = !embree;

                initScene(sceneIndex);
                forceDisplay = true;

                break;
#endif

            case 'w':
                std::cout << "Save config.txt" << std::endl;
                writeConfig("config.txt");
                break;

            case 'o':
                std::cout << "Read config.txt" << std::endl;
                readConfig("config.txt");
                initScene(sceneIndex);
                forceDisplay = true;
                break;

            case 'n':
                normalization = !normalization;

                forceDisplay = true;
                break;

            case 32:
                animated = !animated;
                break;

            case 13:
                drawGeometryType++;
                drawGeometryType = drawGeometryType % 4;
                break;
            }
        }
    private:
        std::vector<float> v0;
        std::vector<float> v1;

        GeometryOccluderSet* occluderSet = nullptr;
        HelperTriangleMeshContainer* meshContainer = nullptr;
        VisualDebugger* debugger = nullptr;
        VisibilityResult result = UNKNOWN;
        int vertexCount = 3;
        bool silhouetteOptimisation = true;
        int drawGeometryType = 0;
        bool representativeLineSampling = true;
        QueryConfiguration::PrecisionType precisionType = QueryConfiguration::DOUBLE;
        bool detectApertureOnly = false;
        bool animated = false;
        bool normalization = true;
        float scaling = 0.1f;
        float phi = 0;
        float eta = 0;
        int   sceneIndex = 2;
        float globalScaling = 1;
        bool forceDisplay = true;
#if EMBREE
        bool embree = false;
#endif
    };
}

static VisilibDemoMain* demo = nullptr;

void display()
{
    demo->display();
    glutSwapBuffers();
}

void keyboard(unsigned char key, int x, int y)
{
    demo->keyboard(key, x, y);
}

void animate()
{
    demo->animate();
}

int main(int argc, char** argv)
{
    demo = new VisilibDemoMain();

    demo->writeHelp();

    glutInit(&argc, argv);

    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    glutInitWindowSize(800, 800);

    glutInitWindowPosition(100, 100);

    glutCreateWindow("Visilib demo");

    glutDisplayFunc(display);

    glutKeyboardFunc(keyboard);

    glutIdleFunc(animate);

    glutReshapeFunc(zprReshape);

    glutMouseFunc(zprMouse);

    glutMotionFunc(zprMotion);

    if (!demo->init())
    {
        std::cout << "Error reading geometry files. Exit" << std::endl;
        return 1;
    }
    glutMainLoop();

    delete demo;

    return 0;
}

#if EMBREE
RTCDevice SilhouetteContainerEmbree::mDevice = nullptr;
#endif
