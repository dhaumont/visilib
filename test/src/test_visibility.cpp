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

#include <string>
#include <iostream>

#include "../demo/demo_helper.h"
#include "test_definition.h"

using namespace visilib;
using namespace visilibDemo;


bool testConfiguration(const DemoConfiguration& configuration, VisibilityResult expected, const std::string& filename)
{
    std::vector<float> v0, v1;


    DemoHelper::generatePolygon(v0, configuration.vertexCount0, configuration.scaling0, configuration.phi - 3.14519f, configuration.globalScaling);
    DemoHelper::generatePolygon(v1, configuration.vertexCount1, configuration.scaling1, configuration.phi, configuration.globalScaling);

    VisibilityExactQueryConfiguration config;
    config.silhouetteOptimization = true;
    config.normalization = configuration.normalization;
    config.precision = configuration.precisionType;
    config.detectApertureOnly = false;

    DemoConfiguration::SCENE_TYPE s = DemoConfiguration::SLOT_OFF_AXIS_01;
    auto meshContainer = DemoHelper::createScene(s, configuration.globalScaling);
    GeometryOccluderSet* occluderSet = DemoHelper::createOccluderSet(meshContainer);

     HelperVisualDebugger debugger;
    debugger.getTimeout().setMaximumDuration(queryBudgetSeconds());
    debugger.setDisplayStatistics(false);

    auto result = visilib::areVisible(occluderSet, &v0[0], v0.size() / 3, &v1[0], v1.size() / 3, config, &debugger);
  
    return (result == expected);
}

bool VisibilityTest(std::string& )
{
    std::vector<VisibilityExactQueryConfiguration::PrecisionType> precisions = { VisibilityExactQueryConfiguration::DOUBLE };
#ifdef ENABLE_LEDA
    precisions.push_back(VisibilityExactQueryConfiguration::LEDA_REAL);
#endif
#ifdef ENABLE_MPFR
    precisions.push_back(VisibilityExactQueryConfiguration::MPFR);
#endif
#ifdef ENABLE_GMP
    precisions.push_back(VisibilityExactQueryConfiguration::GMP_FLOAT);
    precisions.push_back(VisibilityExactQueryConfiguration::GMP_RATIONAL);
#endif
    std::vector<float> globalScalings = {0.05f, 0.1f, 1.f, 5.0f };
    std::vector<size_t> vertexCount = { 1,2,3,5,7,9,11 };
    std::vector<NormalizationType> normalizations = { HYPERSPHERE_NORMALIZATION, NO_NORMALIZATION, CHART_NORMALIZATION };
    std::vector<float> phis = { 0.0};
    std::vector<bool> silhouetteOptimization = {true,false};
    bool result = true;
  
    int i = 0;
    size_t passedCount = 0;
    for (auto silhouette : silhouetteOptimization)
    {
        for (auto precision : precisions)
        {
            for (auto globalScaling : globalScalings)
            {
                for (auto v0 : vertexCount)
                {
                    for (auto v1 : vertexCount)
                    {
                        for (auto phi : phis)
                        {
                            for (auto normalization : normalizations)
                            {
                                DemoConfiguration configuration;
                                configuration.globalScaling = globalScaling;
                                configuration.phi = phi;
                                configuration.vertexCount0 = v0;
                                configuration.vertexCount1 = v1;
                                configuration.solverType = VisibilityExactQueryConfiguration::EXACT_APERTURE_FINDER;
                                configuration.exactArithmetic = false;
                                configuration.normalization = normalization;
                                configuration.scaling0 = 0.14;
                                configuration.scaling1 = 0.14;
                                configuration.silhouetteOptimisation = silhouette;
                                std::stringstream ss;
                                ss << v0 << "_" << v1 << ".obj";

                                if (!testConfiguration(configuration, VISIBLE, ss.str()))
                                {   configuration.displaySummary();
                                    std::cout << " FAILED" << std::endl;
                                    result = false;
                                }
                                else
                                {
                                    passedCount++;
                                }
                                i++;
                            }
                        }
                    }
                }
            }
        }
    }

    std::cout << "VisibilityTest: " << passedCount << " of " << i << " configurations passed" << std::endl;

    return result;
}

namespace
{
    /* A wall standing across the light: the quad of the plane x, between y0 and y1,
    between z0 and z1, as the two triangles of a mesh. */
    HelperTriangleMesh* createWall(float x, float y0, float y1, float z0, float z1)
    {
        std::vector<float> vertices =
        {
            x, y0, z0,
            x, y1, z0,
            x, y1, z1,
            x, y0, z1
        };
        std::vector<int> triangles = { 0, 1, 2, 0, 2, 3 };
        return new HelperTriangleMesh(vertices, triangles);
    }

    GeometryOccluderSet* createWallScene(const std::vector<HelperTriangleMesh*>& walls)
    {
        HelperTriangleMeshContainer* container = new HelperTriangleMeshContainer();
        for (auto wall : walls)
        {
            container->add(wall);
        }
        return DemoHelper::createOccluderSet(container);
    }

    /* One query of the directional pipeline: two unit squares facing each other along
    x, their centers giving the direction of the light, solved by one of the two
    exact solvers. */
    VisibilityResult runDirectionalQuery(GeometryOccluderSet* scene, VisibilityExactQueryConfiguration::SolverType solverType)
    {
        std::vector<float> v0 = { 0, -1, -1, 0, 1, -1, 0, 1, 1, 0, -1, 1 };
        std::vector<float> v1 = { 2, -1, -1, 2, -1, 1, 2, 1, 1, 2, 1, -1 };

        VisibilityExactQueryConfiguration config;
        config.directionalLight = true;
        config.precision = VisibilityExactQueryConfiguration::DOUBLE;
        config.solverType = solverType;
        config.detectApertureOnly = false;

        HelperVisualDebugger debugger;
        debugger.getTimeout().setMaximumDuration(queryBudgetSeconds());
        debugger.setDisplayStatistics(false);

        return visilib::areVisible(scene, &v0[0], v0.size() / 3, &v1[0], v1.size() / 3, config, &debugger);
    }

    bool checkDirectionalCase(const char* name, const std::vector<HelperTriangleMesh*>& walls, VisibilityResult expected)
    {
        bool result = true;
        const VisibilityExactQueryConfiguration::SolverType solvers[] =
        {
            VisibilityExactQueryConfiguration::EXACT_SEQUENTIAL_SOLVER,
            VisibilityExactQueryConfiguration::EXACT_APERTURE_FINDER
        };

        for (auto solverType : solvers)
        {
            GeometryOccluderSet* scene = createWallScene(walls);
            const VisibilityResult answer = runDirectionalQuery(scene, solverType);

            if (answer != expected)
            {
                std::cout << "DirectionalVisibilityTest: " << name << " with solver "
                          << (solverType == VisibilityExactQueryConfiguration::EXACT_SEQUENTIAL_SOLVER ? "SEQUENTIAL" : "APERTURE")
                          << " answered " << (int)answer << ", expected " << (int)expected << std::endl;
                result = false;
            }
        }
        return result;
    }
}

/* The directional pipeline, end to end: the query of VisibilityExactQuery_ run with
the family as its point type, from the factory of visilib.hpp to the answer. The
scenes are the four shapes a wall can take across the light between two squares. */
bool DirectionalVisibilityTest(std::string& )
{
    bool result = true;

    // No occluder: the two squares see each other.
    if (!checkDirectionalCase("empty scene", {}, VISIBLE))
        result = false;

    // A wall wider than both shadows: nothing of the light gets through.
    if (!checkDirectionalCase("full wall", { createWall(1.0f, -3.0f, 3.0f, -3.0f, 3.0f) }, HIDDEN))
        result = false;

    // A wall covering half of the beam: the other half passes.
    if (!checkDirectionalCase("half wall", { createWall(1.0f, -3.0f, 0.0f, -3.0f, 3.0f) }, VISIBLE))
        result = false;

    // A patch outside the shadow of either square: it hides nothing.
    if (!checkDirectionalCase("off beam patch", { createWall(1.0f, 2.0f, 3.0f, 2.0f, 3.0f) }, VISIBLE))
        result = false;

    if (!result)
        return false;

    std::cout << "DirectionalVisibilityTest SUCCESS" << std::endl;
    return true;
}
