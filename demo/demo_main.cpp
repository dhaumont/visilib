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

#include <imgui.h>
#include <backends/imgui_impl_glut.h>
#include <backends/imgui_impl_opengl2.h>

#include <string>
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>

// #include "xmmintrin.h"
// #include "pmmintrin.h"

#include "helper_triangle_mesh_container.h"
#include "silhouette_container_embree.h"
#include "helper_synthetic_mesh_builder.h"
#include "helper_geometry_scene_reader.h"
#include "demo_viewer_glut.h"
#include "demo_debug_visualisation_gl.h"
#include "demo_helper.h"
#include "math_gmp_float.h"
#include "math_mpfr.h"
#include "math_arithmetic.h"
using namespace visilib;
using namespace std;
using namespace visilibDemo;

namespace visilibDemo
{
    class VisilibDemoMain;

    class VisilibDemoMain
    {
    public:
       bool show_demo_window;
       int wc = 10;
       int rc = 0;
        
        VisilibDemoMain()
        {
            show_demo_window = false;
        }
        ~VisilibDemoMain()
        {
        }
        bool init()
        {
//            _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
//            _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
            readConfig("config.txt");
            forceDisplay = true;
#ifdef USE_GLUT
            gluLookAt(1, 1, 1,
                0, 0, 0,
                1, 0, 0);
#endif
            return initScene(mDemoConfiguration.sceneIndex);
        }
        bool initScene(DemoConfiguration::SCENE_TYPE s)
        {
            srand(0);

            delete debugger;
            debugger = new HelperVisualDebugger();

            delete meshContainer;
            meshContainer = DemoHelper::createScene(s, mDemoConfiguration.globalScaling);

            delete occluderSet;
            occluderSet = DemoHelper::createOccluderSet(meshContainer);

            return true;
        }

        void resolveVisibility()
        {
            VisibilityExactQueryConfiguration config;
            config.silhouetteOptimization = mDemoConfiguration.silhouetteOptimisation;
            config.hyperSphereNormalization = mDemoConfiguration.normalization;
            config.precision = mDemoConfiguration.precisionType;
            config.solverType = mDemoConfiguration.solverType;
            config.detectApertureOnly = mDemoConfiguration.detectApertureOnly;
            config.tolerance = mDemoConfiguration.tolerance;
            config.confidenceValue = mDemoConfiguration.confidenceValue;
            double areaSources = 2 * mDemoConfiguration.scaling;
            areaSources *= areaSources;
            double minimumAreAperture =  mDemoConfiguration.minimumApertureSize;
            minimumAreAperture *= minimumAreAperture;
            config.minimumNormalizedApertureSize = minimumAreAperture/ areaSources;
#if EMBREE 
            config.useEmbree = mDemoConfiguration.embree;
#endif
#if ENABLE_MPFR
            MathMpfr::set_default_precision(MathArithmetic<MathMpfr>::bitsCount());
#endif
#if ENABLE_GMP
            MathGmpFloat::set_default_precision(MathArithmetic<MathGmpFloat>::bitsCount());
#endif


            result = visilib::areVisible(occluderSet, &v0[0], v0.size() / 3, &v1[0], v1.size() / 3, config, debugger);
        }

        void animate()
        {
            if (animated)
            {
                mDemoConfiguration.phi += 0.005f;
                mDemoConfiguration.eta += 0.001f;
                forceDisplay = true;
            }

            if (forceDisplay)
            {
                DemoHelper::generatePolygon(v0, mDemoConfiguration.vertexCount0, mDemoConfiguration.scaling, mDemoConfiguration.phi - (float)M_PI, mDemoConfiguration.globalScaling);
                DemoHelper::generatePolygon(v1, mDemoConfiguration.vertexCount1, mDemoConfiguration.scaling, mDemoConfiguration.phi, mDemoConfiguration.globalScaling);

                resolveVisibility();

                forceDisplay = false;
            }
#ifdef USE_GLUT
            glutPostRedisplay();
#endif
        }

        void writeConfig(const std::string & filename)
        {
            mDemoConfiguration.writeConfig(filename);
        }

        void readConfig(const std::string & filename)
        {
            mDemoConfiguration.readConfig(filename);
        }
        
        void showGUI()
        {
           if (show_demo_window)
            {
                ImGui::ShowDemoWindow(&show_demo_window);
            }

            if (ImGui::BeginMainMenuBar()) {
                if (ImGui::BeginMenu("File")) {
                    if (ImGui::MenuItem("Open...", NULL, false, true)) {
                        readConfig("config.txt");
                        forceDisplay = true;
                    }
                    if (ImGui::MenuItem("Save...", NULL, false, true)) {
                        writeConfig("config.txt");
                    }
                    if (ImGui::MenuItem("Quit", NULL, false, true)) {
                        exit(0);
                    }
                  
                    ImGui::EndMenu();
                }
                ImGui::EndMainMenuBar();
            }

           

            {
                ImGui::Begin("Visilib 1.0. Demo application");

                forceDisplay |= ImGui::Checkbox("silhouette optimisation", &mDemoConfiguration.silhouetteOptimisation);
                forceDisplay |= ImGui::Checkbox("normalization", &mDemoConfiguration.normalization);

                const char* resolve_mode_items[] = { "Recursive", "NonRecursive", "Compare" };
               
             
#if EMBREE
                if (ImGui::Checkbox("embree ray tracing", mDemoConfiguration.embree)) {
                    initScene(mDemoConfiguration.sceneIndex);
                    forceDisplay = true;
                }
#endif
                //forceDisplay |= ImGui::Checkbox("representative line sampling strategy", &mDemoConfiguration.representativeLineSampling);
                forceDisplay |= ImGui::Checkbox("detect aperture only", &mDemoConfiguration.detectApertureOnly);

            
                const char* precision_items[] = { 
                    "FLOAT",     
                    "DOUBLE"  
            #ifdef ENABLE_LEDA
                    , "LEDA_REAL" 
        #endif
        #ifdef ENABLE_GMP
                    , "GMP_FLOAT"    
                    , "GMP_RATIONAL"
        #endif
        #ifdef ENABLE_MPFR
                    , "MPFR"
        #endif
                };
                int precisionType = mDemoConfiguration.precisionType;
                forceDisplay |= ImGui::Combo("Precision", &precisionType, precision_items, IM_ARRAYSIZE(precision_items));
                mDemoConfiguration.precisionType = (VisibilityExactQueryConfiguration::PrecisionType)precisionType;

                ImGui::SameLine();
                ImGui::Text(DemoConfiguration::toStr(mDemoConfiguration.precisionType).c_str()); // validate that the enum and the string list match up...

            
                ImGui::SameLine();
                ImGui::Text("%i", mDemoConfiguration.sceneIndex);
                ImGui::SameLine();
                if (ImGui::Button("next ")) {
                    if (mDemoConfiguration.sceneIndex < 9) 
                    {
                  //      mDemoConfiguration.sceneIndex ++;
                    }
                    else
                    {
                  //      mDemoConfiguration.sceneIndex = 0;
                    }
                    initScene(mDemoConfiguration.sceneIndex);
                    forceDisplay = true;
                }

                forceDisplay |= ImGui::RadioButton("Edit S",&active_source,0);
                forceDisplay |= ImGui::RadioButton("Edit R",&active_source,1);
                forceDisplay |= ImGui::RadioButton("Edit R & S",&active_source,2);
             
                int vertexCount0 = mDemoConfiguration.vertexCount0;
                int vertexCount1 = mDemoConfiguration.vertexCount1;
                int vertexCount2 = mDemoConfiguration.vertexCount1;
           
                float scaling0 = mDemoConfiguration.scaling;
                float scaling1 = mDemoConfiguration.scaling;
                float scaling2 = mDemoConfiguration.scaling;
           
                if (active_source == 0 )
                {
                    forceDisplay |= ImGui::SliderInt("Vertices S", &vertexCount0, 1, 12);
                    
                    forceDisplay |= ImGui::SliderFloat("Size S", &scaling0, 0.01, 1.0);
                    mDemoConfiguration.vertexCount0 = vertexCount0;
                    mDemoConfiguration.scaling = scaling0;
                }
                else if (active_source == 1)
                {
                    forceDisplay |= ImGui::SliderInt("Vertices R", &vertexCount1, 1, 12);
                    forceDisplay |= ImGui::SliderFloat("Size R", &scaling1, 0.01, 1.0);
                    
                     mDemoConfiguration.vertexCount1 = vertexCount1;
                      mDemoConfiguration.scaling = scaling1;
                }
                else if (active_source == 2)
                {
                    forceDisplay |= ImGui::SliderInt("Vertices R & S", &vertexCount2, 1, 12);
                    forceDisplay |= ImGui::SliderFloat("Size S", &scaling2, 0.01, 1.0);
                    
                    mDemoConfiguration.vertexCount0 = vertexCount2;
                    mDemoConfiguration.vertexCount1 = vertexCount2;
                     mDemoConfiguration.scaling = scaling2;
                }
             

                ImGui::Text("adjust global scaling %i%%", (int)round(mDemoConfiguration.globalScaling * 100.0f)); ImGui::SameLine();
                if (ImGui::Button("+")) {
                    mDemoConfiguration.globalScaling *= 2;
                    forceDisplay = true;
                    setViewPortScaling(mDemoConfiguration.globalScaling);
                    initScene(mDemoConfiguration.sceneIndex);
                }
                ImGui::SameLine();
                if (ImGui::Button("-")) {
                    mDemoConfiguration.globalScaling /= 2;
                    forceDisplay = true;
                    setViewPortScaling(mDemoConfiguration.globalScaling);
                    initScene(mDemoConfiguration.sceneIndex);
                }

                if (ImGui::Button("change geometry mode")) {
                    drawGeometryType++;
                    drawGeometryType = drawGeometryType % 4;
                    forceDisplay = true;
                }
                ImGui::SameLine();
                ImGui::Text("%i", drawGeometryType);

                if (ImGui::Button("step")) {
                //    step = true;
                    animated = false;
                    forceDisplay = true;
                }
                ImGui::SameLine();
                forceDisplay |= ImGui::Checkbox("animation", &animated);

                if (ImGui::Button("write config")) {
                    std::stringstream ss;

                    std::cout << "Save config.txt" << std::endl;
                    ss << "config_" << wc++ << ".txt";
                    writeConfig(ss.str());
                }
                ImGui::SameLine();
                if (ImGui::Button("open config")) {
                    std::stringstream ss;

                    std::cout << "Read config.txt" << std::endl;
                    ss << "config_" << rc++ << ".txt";
                    if (rc > wc)
                        rc = 0;
                    readConfig(ss.str());
                    initScene(mDemoConfiguration.sceneIndex);
                }

                ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
                ImGui::End();
            }

            {
                ImGui::Begin("Output");
               // ImGui::Text(output.c_str());
                ImGui::End();
            }

            { 
                std::map<std::string, double> timings = debugger->getTimings();
                std::map<std::string, int> counters = debugger->getCounters();
                ImGui::Begin("Stats");                
                for (const auto& timing : timings)
                {
                    ImGui::Text("%s : %.2f%%", timing.first.c_str(), timing.second);
                }
                for (const auto& counter : counters)
                { 
                    ImGui::Text("%s : %i", counter.first.c_str(), counter.second);
                }
           
                ImGui::End();
            }
        }

        void display()
        {
            DemoDebugVisualisationGl::display(debugger, *meshContainer, v0, v1, result, drawGeometryType,mDemoConfiguration.minimumApertureSize * mDemoConfiguration.globalScaling);
        }

        void displaySettings()
        {
            mDemoConfiguration.displaySettings();
        }

        void writeHelp()
        {
            std::cout << "Visilib 1.0. Demo application " << std::endl;

            std::cout << "  s: enable/disable silhouette optimisation" << std::endl;
            std::cout << "  n: enable/disable nomalization" << std::endl;
            std::cout << "  e: enable/disable exact arithmetic" << std::endl;

#if EMBREE
            std::cout << "  g: enable/disable embree ray tracing" << std::endl;
#endif
            std::cout << "  f: enable/disable detect aperture only" << std::endl;

            std::cout << "  a: enable/disable aggressive" << std::endl;

            std::cout << "  n: enable/disable Plucker normalization" << std::endl;
            
            std::cout << "  x: change scene " << std::endl;
            std::cout << "  +/-: increase/decrease scaling of query polygons" << std::endl;
            std::cout << "  1/2: increase/decrease number of vertices of query polygons" << std::endl;
            std::cout << "  /-*: increase/decrease global scaling" << std::endl;
            std::cout << "  [/]: increase/decrease scaling" << std::endl;

            std::cout << "  w: write config" << std::endl;
            std::cout << "  o: open config" << std::endl;

            std::cout << "  space: start/pause animation" << std::endl;
            std::cout << "  Enter: show/hide geometry" << std::endl;

            std::cout << "  h: write this help" << std::endl;
        }

        void keyboard(unsigned char key, int, int)
        {           
            std::stringstream ss;
            int current = 0;
            switch (key)
            {
            case 27:  // The escape key
            case 'Q':
            case 'q':
                exit(0);   // Simply exit
                break;

            case '>':
                if (mDemoConfiguration.vertexCount1 < 12)
                    mDemoConfiguration.vertexCount1++;

                forceDisplay = true;
                break;
            case '<':
                if (mDemoConfiguration.vertexCount1 > 1)
                    mDemoConfiguration.vertexCount1--;

                forceDisplay = true;
                break;

            case '+':
                if (mDemoConfiguration.scaling < 1.00f)
                    mDemoConfiguration.scaling += 0.01f;
                forceDisplay = true;
                break;

            case '-':
                if (mDemoConfiguration.scaling > 0.02f)
                    mDemoConfiguration.scaling -= 0.01f;

                forceDisplay = true;
                break;

            case '*':
                mDemoConfiguration.globalScaling *= 2;
                forceDisplay = true;
                setViewPortScaling(mDemoConfiguration.globalScaling);
                initScene(mDemoConfiguration.sceneIndex);

                break;

            case '/':
                mDemoConfiguration.globalScaling /= 2;

                forceDisplay = true;
                setViewPortScaling(mDemoConfiguration.globalScaling);
                initScene(mDemoConfiguration.sceneIndex);

                break;

            case ']':                
               if (mDemoConfiguration.minimumApertureSize < 0.1f)
                    mDemoConfiguration.minimumApertureSize *= 1.3;
                forceDisplay = true;                                

                break;

            case '[':
                if (mDemoConfiguration.minimumApertureSize > 0.00001f)
                    mDemoConfiguration.minimumApertureSize /= 1.3;
                forceDisplay = true;                
                break;

            case 'h':
                writeHelp();

                displaySettings();
                break;
            case 's':
                mDemoConfiguration.silhouetteOptimisation = !mDemoConfiguration.silhouetteOptimisation;
                forceDisplay = true;
                break;
            case 'a':
                current = static_cast<int>(mDemoConfiguration.solverType);
                current++;
                if (current >= static_cast<int>(VisibilityExactQueryConfiguration::SOLVER_COUNT)) current = VisibilityExactQueryConfiguration::EXACT_APERTURE_FINDER;
                mDemoConfiguration.solverType = static_cast<VisibilityExactQueryConfiguration::SolverType>(current);
                displaySettings();
                forceDisplay = true;
               break;
            case 'f':
                mDemoConfiguration.detectApertureOnly = !mDemoConfiguration.detectApertureOnly;
                forceDisplay = true;

                break;

            case '0': mDemoConfiguration.sceneIndex = DemoConfiguration::SLOT_OFF_AXIS_01; forceDisplay = true; initScene(mDemoConfiguration.sceneIndex); break;
            case '1': mDemoConfiguration.sceneIndex = DemoConfiguration::SLOT_OFF_AXIS_02; forceDisplay = true;  initScene(mDemoConfiguration.sceneIndex);break;
            case '2': mDemoConfiguration.sceneIndex = DemoConfiguration::SLOT_OFF_AXIS_03; forceDisplay = true;  initScene(mDemoConfiguration.sceneIndex);break;
            case '3': mDemoConfiguration.sceneIndex = DemoConfiguration::SLOT_OFF_AXIS_04; forceDisplay = true; initScene(mDemoConfiguration.sceneIndex); break;
            case '4': mDemoConfiguration.sceneIndex = DemoConfiguration::SLOT_OFF_AXIS_05; forceDisplay = true; initScene(mDemoConfiguration.sceneIndex); break;
            case '5': mDemoConfiguration.sceneIndex = DemoConfiguration::SLOT_ON_AXIS_01; forceDisplay = true;  initScene(mDemoConfiguration.sceneIndex);break;
            case '6': mDemoConfiguration.sceneIndex = DemoConfiguration::SPHERE_WITH_HOLES; forceDisplay = true; initScene(mDemoConfiguration.sceneIndex); break;
            case '7': mDemoConfiguration.sceneIndex = DemoConfiguration::NOISY_SPHERE; forceDisplay = true; initScene(mDemoConfiguration.sceneIndex); break;
            case '8': mDemoConfiguration.sceneIndex = DemoConfiguration::SIMPLE_CUBE; forceDisplay = true; initScene(mDemoConfiguration.sceneIndex); break;
            case '9': mDemoConfiguration.sceneIndex = DemoConfiguration::NOISY_CUBE; forceDisplay = true; initScene(mDemoConfiguration.sceneIndex); break;

            case 'e':
                current = static_cast<int>(mDemoConfiguration.precisionType);
                current++;
                if (current >= static_cast<int>(VisibilityExactQueryConfiguration::PRECISION_COUNT)) current = VisibilityExactQueryConfiguration::FLOAT;
                mDemoConfiguration.precisionType = static_cast<VisibilityExactQueryConfiguration::PrecisionType>(current);
                displaySettings();
                forceDisplay = true;
                break;
#if EMBREE
            case 'g':
                mDemoConfiguration.embree = !mDemoConfiguration.embree;

                initScene(mDemoConfiguration.sceneIndex);
                forceDisplay = true;

                break;
#endif

            case 'w':
                std::cout << "Save config.txt" << std::endl;
                ss << "config_" << wc++ << ".txt";
                writeConfig(ss.str());
                break;

            case 'o':
                std::cout << "Read config.txt" << std::endl;
                ss << "config_" << rc++ << ".txt";
                if (rc > wc)
                    rc = 0;
                readConfig(ss.str());
                initScene(mDemoConfiguration.sceneIndex);
                forceDisplay = true;
                break;

            case 'n':
                mDemoConfiguration.normalization = !mDemoConfiguration.normalization;

                forceDisplay = true;
                break;

            case 32:
                show_demo_window = !show_demo_window;
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

        DemoConfiguration mDemoConfiguration;
        bool forceDisplay = true;
        bool animated = false;
        int drawGeometryType = 0;
        int active_source = 0;
    };
}

static VisilibDemoMain* demo = nullptr;

void display()
{
    ImGui_ImplOpenGL2_NewFrame();
    ImGui_ImplGLUT_NewFrame();
    demo->display();

    ImGui::NewFrame();
    ImGuiIO& io = ImGui::GetIO();
    demo->showGUI();

    ImGui::Render();
    glViewport(0, 0, (GLsizei)io.DisplaySize.x, (GLsizei)io.DisplaySize.y);
    ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
#ifdef USE_GLUT
    glutSwapBuffers();
#endif
}

void keyboard(unsigned char key, int x, int y)
{
#if 1
    ImGui_ImplGLUT_KeyboardFunc(key, x, y);

    ImGuiIO& io = ImGui::GetIO();
    if (!io.WantCaptureKeyboard)
#endif
    {
        demo->keyboard(key, x, y);
    }
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

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls

    // Setup Dear ImGui style
    //ImGui::StyleColorsDark();
    ImGui::StyleColorsClassic();

    // Setup Platform/Renderer bindings
    ImGui_ImplGLUT_Init();
    ImGui_ImplGLUT_InstallFuncs();
    ImGui_ImplOpenGL2_Init();

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
#if 1
    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGLUT_Shutdown();
    ImGui::DestroyContext();
#endif
    return 0;
}

