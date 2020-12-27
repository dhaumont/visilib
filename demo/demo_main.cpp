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
            return initScene(configuration.sceneIndex);
        }
        bool initScene(int s)
        {
            srand(0);

            delete debugger;
            debugger = new HelperVisualDebugger();

            delete meshContainer;
            meshContainer = DemoHelper::createScene(s, configuration.globalScaling);
            
            delete occluderSet;
            occluderSet = DemoHelper::createOccluderSet(meshContainer);

            return true;
        }

        void resolveVisibility()
        {
            VisibilityExactQueryConfiguration config;
            config.silhouetteOptimization = configuration.silhouetteOptimisation;
            config.hyperSphereNormalization = configuration.normalization;
            config.precision = configuration.precisionType;
            config.representativeLineSampling = configuration.representativeLineSampling;
            config.detectApertureOnly = configuration.detectApertureOnly;
#if EMBREE 
            config.useEmbree = embree;
#endif
            result = visilib::areVisible(occluderSet, &v0[0], v0.size() / 3, &v1[0], v1.size() / 3, config, debugger);
        }

        void animate()
        {
            if (animated)
            {
                configuration.phi += 0.005f;
                configuration.eta += 0.001f;
                forceDisplay = true;
            }

            if (forceDisplay)
            {
                DemoHelper::generatePolygon(v0, configuration.vertexCount0, configuration.scaling, configuration.phi - 3.14519f, configuration.globalScaling);
                DemoHelper::generatePolygon(v1, configuration.vertexCount1, configuration.scaling, configuration.phi, configuration.globalScaling);

                resolveVisibility();

                forceDisplay = false;
            }
#ifdef USE_GLUT
            glutPostRedisplay();
#endif
        }

        void writeConfig(const std::string & filename)
        {
            configuration.writeConfig(filename);
        }

        void readConfig(const std::string & filename)
        {
            configuration.readConfig(filename);
        }

        void display()
        {
            DemoDebugVisualisationGl::display(debugger, *meshContainer, v0, v1, result, drawGeometryType);
        }

        void displaySettings()
        {
            configuration.displaySettings();
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
                if (configuration.vertexCount1 < 12)
                    configuration.vertexCount1++;

                forceDisplay = true;
                break;
            case '1':
                if (configuration.vertexCount1 > 1)
                    configuration.vertexCount1--;

                forceDisplay = true;
                break;

            case '+':
                if (configuration.scaling < 1.00f)
                    configuration.scaling += 0.01f;
                forceDisplay = true;
                break;

            case '-':
                if (configuration.scaling > 0.02f)
                    configuration.scaling -= 0.01f;

                forceDisplay = true;
                break;

            case '*':
                configuration.globalScaling *= 2;
                forceDisplay = true;
                setViewPortScaling(configuration.globalScaling);
                initScene(configuration.sceneIndex);

                break;

            case '/':
                configuration.globalScaling /= 2;

                forceDisplay = true;
                setViewPortScaling(configuration.globalScaling);
                initScene(configuration.sceneIndex);

                break;

            case 'h':

                writeHelp();

                displaySettings();
                break;
            case 's':
                configuration.silhouetteOptimisation = !configuration.silhouetteOptimisation;
                forceDisplay = true;
                break;
            case 'f':
                configuration.detectApertureOnly = !configuration.detectApertureOnly;
                forceDisplay = true;

                break;

            case 'x':
                configuration.sceneIndex++;
                if (configuration.sceneIndex > 9)
                    configuration.sceneIndex = 0;
                initScene(configuration.sceneIndex);
                forceDisplay = true;
                break;
            case 'r':
                configuration.representativeLineSampling = !configuration.representativeLineSampling;
                forceDisplay = true;
                break;


#ifdef EXACT_ARITHMETIC
            case 'e':
                precisionType = precisionType == VisibilityExactQueryConfiguration::DOUBLE ? VisibilityExactQueryConfiguration::EXACT : VisibilityExactQueryConfiguration::DOUBLE;
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
                initScene(configuration.sceneIndex);
                forceDisplay = true;
                break;

            case 'n':
                configuration.normalization = !configuration.normalization;

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

        DemoConfiguration configuration;
        bool forceDisplay = true;
        bool animated = false;
        int drawGeometryType = 0;


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
