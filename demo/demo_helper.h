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

#pragma once

#include <vector>
#include "visilib.h"
#include "helper_triangle_mesh_container.h"

namespace visilibDemo
{
    struct DemoConfiguration
    {
        enum SCENE_TYPE
        {
            SINGLE_TRIANGLE_01,
            SINGLE_TRIANGLE_02,
            SINGLE_TRIANGLE_03,
            SINGLE_TRIANGLE_04,
            SLOT_OFF_AXIS_01,
            SLOT_OFF_AXIS_02,
            SLOT_OFF_AXIS_03,
            SLOT_OFF_AXIS_04,
            SLOT_OFF_AXIS_05,
            SLOT_ON_AXIS_01,
            SPHERE_WITH_HOLES,
            REGULAR_GRID,
            NOISY_SPHERE,
            SIMPLE_CUBE,
            NOISY_CUBE,
            CUBES_AND_GRID,
            CUBES_AND_GRID_100,
            COUNT
        };

        std::string getStatusString(bool enable)
        {
            return enable ? "ON" : "OFF";
        }
        size_t vertexCount0 = 3;
        size_t vertexCount1 = 3;
        bool silhouetteOptimisation = true;
        VisibilityExactQueryConfiguration::PrecisionType precisionType = VisibilityExactQueryConfiguration::DOUBLE;
        bool detectApertureOnly = false;
        bool normalization = false;
        float scaling = 0.1f;
        float phi = 0;
        float eta = 0;
        DemoConfiguration::SCENE_TYPE sceneIndex = DemoConfiguration::SINGLE_TRIANGLE_01;
        float globalScaling = 1;
        double tolerance = -1;
        double minimumApertureSize = 0.01;
        double confidenceValue = 0.99;
        VisibilityExactQueryConfiguration::SolverType solverType = VisibilityExactQueryConfiguration::EXACT_APERTURE_FINDER;
        bool exactArithmetic = false;

        /* What the solver consistency test expects the solvers to answer, "VISIBLE" or
        "HIDDEN", empty when the case carries no expectation. Test data living with the
        configuration it pins: the result of a hand check travels with the case file into
        the test, instead of a table that has to be kept in sync with the folder. */
        std::string expected;
#if EMBREE
        bool embree = false;
#endif
        void displaySettings();
        void readConfig(const std::string& filename);
        void writeConfig(const std::string& filename);

        /* The settings as a JSON document, without touching the disk. */
        std::string jsonString() const;

        /* The name these settings give themselves: the scene and every parameter that
        changes what is computed, so that two different setups never land on the same
        file. Used for the test case files. */
        std::string caseFileName() const;

        /* Stores the settings in the test case folder under caseFileName(). Never
        overwrites: a file of that name already holds these settings, or it was edited
        since and must not be destroyed. The reason ends up in message. */
        bool writeTestCase(const std::string& folder, std::string& message) const;

        static const char* sceneName(SCENE_TYPE scene);

        void displaySummary();
        static const std::string toStr(VisibilityExactQueryConfiguration::PrecisionType p);
        static const std::string toStr(VisibilityExactQueryConfiguration::SolverType p);
        static const char* resultName(visilib::VisibilityResult result);
        static visilib::VisibilityResult resultFromName(const std::string& name);
    };

    class DemoHelper
    {
    public:
        static void generatePolygon(std::vector<float>& v, size_t vertexCount, float size, float phi, float scaling);
        static bool load(visilib::HelperTriangleMeshContainer* scene, const std::string& fileName, bool removeDegeneratedTriangles = false);
        static double getApertureSize(DemoConfiguration::SCENE_TYPE s);
        static int getSubdivision(DemoConfiguration::SCENE_TYPE s);
        static visilib::HelperTriangleMeshContainer* createScene(DemoConfiguration::SCENE_TYPE s, float globalScalingFactor);
        static visilib::GeometryOccluderSet* createOccluderSet(visilib::HelperTriangleMeshContainer* aContainer);
        static void configureDemoConfiguration(const std::string& name, DemoConfiguration& configuration);
        static void exportQueryToObj(const std::string& filename, const std::vector<float>& v0, const std::vector<float>& v1, const visilib::HelperTriangleMeshContainer& aContainer, float aScaling = 1.);
        static void appendPolygonRepresentationToFileObj(std::ofstream& stream, int& anOffset, const std::vector<float>& v, float aScaling = 1);

        /* The folder the test cases live in, handed to the demo at configure time, and
        the .json files it holds, sorted by name. Listing a missing folder is not an
        error: it simply has no cases yet. */
        static std::string testCaseFolder();
        static std::vector<std::string> listTestCaseFiles(const std::string& folder);

        /* The short code a precision goes by in file names and result lines. */
        static std::string precisionTag(VisibilityExactQueryConfiguration::PrecisionType precision);
    };
}

