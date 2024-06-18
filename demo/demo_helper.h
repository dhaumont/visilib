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
        bool detectApertureOnly = false;
        bool normalization = true;
        float scaling = 0.1f;
        float phi = 0;
        float eta = 0;
        DemoConfiguration::SCENE_TYPE sceneIndex = DemoConfiguration::SPHERE_WITH_HOLES;
        float globalScaling = 1;
        double tolerance = -1;
        double minimumApertureSize = 0.00175;
        double confidenceValue = 0.99;
        bool sampling = false;
        bool exactArithmetic = false;       
#if EMBREE
        bool embree = false;
#endif
        void displaySettings();
        void readConfig(const std::string& filename);
        void writeConfig(const std::string& filename);
        void displaySummary();        
        visilib::VisibilityExactQueryConfiguration::PrecisionType getPrecisionType() const;        
    };

    class DemoHelper
    {
    public:
        static void generatePolygon(std::vector<float>& v, size_t vertexCount, float size, float phi, float scaling);
        static bool load(visilib::HelperTriangleMeshContainer* scene, const std::string& fileName, bool removeDegeneratedTriangles = false);
        static double getApertureSize(DemoConfiguration::SCENE_TYPE s);
        static visilib::HelperTriangleMeshContainer* createScene(DemoConfiguration::SCENE_TYPE s, float globalScalingFactor);
        static visilib::GeometryOccluderSet* createOccluderSet(visilib::HelperTriangleMeshContainer* aContainer);
        static void configureDemoConfiguration(const std::string& name, DemoConfiguration& configuration);
        static void exportQueryToObj(const std::string& filename, const std::vector<float>& v0, const std::vector<float>& v1, const visilib::HelperTriangleMeshContainer& aContainer, float aScaling = 1.);
        static void appendPolygonRepresentationToFileObj(std::ofstream& stream, int& anOffset, const std::vector<float>& v, float aScaling = 1);
    };
}

