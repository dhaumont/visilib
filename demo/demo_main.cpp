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
#ifdef USE_GLUT
#include <GL/gl.h>
#include <GL/glut.h>
#endif
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
#ifdef USE_GLUT
            gluLookAt(1, 1, 1,
                0, 0, 0,
                1, 0, 0);
#endif
            return initScene(sceneIndex);
        }
        bool initScene(int s)
        {
            srand(0);

            delete debugger;
            debugger = new HelperVisualDebugger();

            delete meshContainer;
            meshContainer = DemoHelper::createScene(s, globalScaling);
            
            delete occluderSet;
            occluderSet = createOccluderSet();

            return true;
        }

        GeometryOccluderSet* createOccluderSet()
        {
            GeometryOccluderSet* occluderSet = new GeometryOccluderSet();
            for (size_t index = 0; index < meshContainer->getGeometryCount(); index++)
            {
                GeometryDiscreteMeshDescription* info = meshContainer->createTriangleMeshDescription(index);
                occluderSet->addOccluder(info);
            }
            occluderSet->prepare();

            return occluderSet;
        }

        void resolveVisibility()
        {
            VisibilityExactQueryConfiguration config;
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
#ifdef USE_GLUT
            glutPostRedisplay();
#endif
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
                else if (tokens[0] == "precisionType") { precisionType = (VisibilityExactQueryConfiguration::PrecisionType)atoi(tokens[2].c_str()); }

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
        HelperVisualDebugger* debugger = nullptr;
        VisibilityResult result = UNKNOWN;
        int vertexCount = 3;
        bool silhouetteOptimisation = true;
        int drawGeometryType = 0;
        bool representativeLineSampling = true;
        VisibilityExactQueryConfiguration::PrecisionType precisionType = VisibilityExactQueryConfiguration::DOUBLE;
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
#ifdef USE_GLUT
    glutSwapBuffers();
#endif
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

#ifdef USE_GLUT
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
#endif
    if (!demo->init())
    {
        std::cout << "Error reading geometry files. Exit" << std::endl;
        return 1;
    }
#ifdef USE_GLUT
    glutMainLoop();
#else
    animate();
#endif
    delete demo;

    return 0;
}

#if EMBREE
RTCDevice SilhouetteContainerEmbree::mDevice = nullptr;
#endif
