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

#include <string>
#include <iostream>

#include "../demo/demo_helper.h"
#include "test_definition.h"

using namespace visilib;
using namespace visilibDemo;


bool testConfiguration(const DemoConfiguration& configuration, VisibilityResult expected
    )
{
    std::vector<float> v0, v1;
    
   
    DemoHelper::generatePolygon(v0, configuration.vertexCount0, configuration.scaling, configuration.phi - 3.14519f, configuration.globalScaling);
    DemoHelper::generatePolygon(v1, configuration.vertexCount1, configuration.scaling, configuration.phi, configuration.globalScaling);

    VisibilityExactQueryConfiguration config;
    config.silhouetteOptimization = true;
    config.hyperSphereNormalization = configuration.normalization;
    config.precision = configuration.precisionType;
    config.representativeLineSampling = true;
    config.detectApertureOnly = false;

    int s = 2;
    auto meshContainer = DemoHelper::createScene(s, configuration.globalScaling);
    GeometryOccluderSet* occluderSet = DemoHelper::createOccluderSet(meshContainer);

    auto result = visilib::areVisible(occluderSet, &v0[0], v0.size() / 3, &v1[0], v1.size() / 3, config, nullptr);
    return (result == expected);
}

bool VisibilityTest(std::string& )
{  
    std::vector<VisibilityExactQueryConfiguration::PrecisionType> precisions = { VisibilityExactQueryConfiguration::DOUBLE };
#ifdef EXACT_ARITHMETIC
    precisions.push_back(VisibilityExactQueryConfiguration::EXACT);
#endif
    std::vector<float> globalScalings = {0.05f, 0.1f, 1.f, 10, 100, 100 };
    std::vector<size_t> vertexCount = { 1,2,3,5,7,9,11 };
    std::vector<bool> normalizations = { true,false };
    std::vector<float> phis = { 0.0};
    bool result = true;
  
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
                            configuration.precisionType = precision;
                            configuration.normalization = normalization;
                            configuration.scaling = 0.14;
                            configuration.displaySummary();

                            if (!testConfiguration(configuration, VISIBLE))
                            {
                                std::cout << " FAILED" << std::endl;
                                result = false;
                            }
                            else
                            {
                                std::cout << " SUCCESS" << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
 
    return result;
}
