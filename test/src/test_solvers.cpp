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
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <future>
#include <iterator>
#include <memory>
#include <thread>
#include <utility>

#include "../demo/demo_helper.h"
#include "test_definition.h"

#include "geometry_convex_polygon.h"
#include "geometry_occluder_set.h"
#include "math_arithmetic.h"
#include "math_geometry.h"
#include "silhouette_mesh_face.h"
#include "math_vector_3.h"

using namespace visilib;
using namespace visilibDemo;

namespace
{
    const size_t cResultWidth = 8; 
    const float cSourceHalfTurn = 3.14519f;
    
    VisibilityResult solveWith(GeometryOccluderSet* occluderSet, const std::vector<float>& v0, const std::vector<float>& v1,
                               const DemoConfiguration& configuration, VisibilityExactQueryConfiguration::SolverType solverType)
    {
        VisibilityExactQueryConfiguration query;
        query.detectApertureOnly = configuration.detectApertureOnly;
        query.hyperSphereNormalization = configuration.normalization;
        query.silhouetteOptimization = configuration.silhouetteOptimisation;
        query.precision = configuration.precisionType;
        query.tolerance = configuration.tolerance;
        query.confidenceValue = configuration.confidenceValue;
        query.solverType = solverType;

        const double side = 2.0 * configuration.scaling;
        const double aperture = configuration.minimumApertureSize;
        query.minimumNormalizedApertureSize = (aperture * aperture) / (side * side);

        return visilib::areVisible(occluderSet, &v0[0], v0.size() / 3, &v1[0], v1.size() / 3, query, nullptr);
    }

    struct SweepVariant
    {
        VisibilityExactQueryConfiguration::PrecisionType mPrecision = VisibilityExactQueryConfiguration::DOUBLE;
        float mGlobalScaling = 1.0f;
        bool mSilhouetteOptimisation = false;
        bool mDetectApertureOnly = false;
    };

    struct SweepBounds
    {
        std::vector<VisibilityExactQueryConfiguration::PrecisionType> mPrecisions;
        std::vector<float> mGlobalScalings;
        std::vector<bool> mSilhouetteOptimisations;
        std::vector<bool> mDetectApertureOnly;

        size_t settingsCount() const
        {
            return mGlobalScalings.size() * mSilhouetteOptimisations.size() * mDetectApertureOnly.size();
        }

        size_t variantCount() const
        {
            return settingsCount() * mPrecisions.size();
        }

        struct Iterator
        {
            typedef std::forward_iterator_tag iterator_category;
            typedef SweepVariant value_type;
            typedef std::ptrdiff_t difference_type;
            typedef const SweepVariant* pointer;
            typedef const SweepVariant& reference;

            Iterator() : mSweep(nullptr), mPrecisionIndex(0), mScalingIndex(0), mSilhouetteIndex(0), mEarlyIndex(0)
            {
            }

            Iterator(const SweepBounds* sweep, size_t precisionIndex, size_t scalingIndex, size_t silhouetteIndex, size_t earlyIndex)
                : mSweep(sweep), mPrecisionIndex(precisionIndex), mScalingIndex(scalingIndex), mSilhouetteIndex(silhouetteIndex), mEarlyIndex(earlyIndex)
            {
            }

            SweepVariant at(size_t precisionIndex, size_t scalingIndex, size_t silhouetteIndex, size_t earlyIndex) const
            {
                SweepVariant variant;
                variant.mPrecision = mSweep->mPrecisions[precisionIndex];
                variant.mGlobalScaling = mSweep->mGlobalScalings[scalingIndex];
                variant.mSilhouetteOptimisation = mSweep->mSilhouetteOptimisations[silhouetteIndex];
                variant.mDetectApertureOnly = mSweep->mDetectApertureOnly[earlyIndex];

                return variant;
            }

            reference operator*() const
            {
                mCached = at(mPrecisionIndex, mScalingIndex, mSilhouetteIndex, mEarlyIndex);
                return mCached;
            }

            pointer operator->() const
            {
                mCached = at(mPrecisionIndex, mScalingIndex, mSilhouetteIndex, mEarlyIndex);
                return &mCached;
            }

            Iterator& operator++()
            {
                mEarlyIndex++;

                if (mEarlyIndex == mSweep->mDetectApertureOnly.size())
                {
                    mEarlyIndex = 0;
                    mSilhouetteIndex++;
                }

                if (mSilhouetteIndex == mSweep->mSilhouetteOptimisations.size())
                {
                    mSilhouetteIndex = 0;
                    mScalingIndex++;
                }

                if (mScalingIndex == mSweep->mGlobalScalings.size())
                {
                    mScalingIndex = 0;
                    mPrecisionIndex++;
                }

                return *this;
            }

            Iterator operator++(int)
            {
                Iterator copy = *this;
                ++*this;
                return copy;
            }

            bool operator==(const Iterator& other) const
            {
                return mSweep == other.mSweep && mPrecisionIndex == other.mPrecisionIndex && mScalingIndex == other.mScalingIndex &&
                       mSilhouetteIndex == other.mSilhouetteIndex && mEarlyIndex == other.mEarlyIndex;
            }

            bool operator!=(const Iterator& other) const
            {
                return !(*this == other);
            }

            const SweepBounds* mSweep;
            size_t mPrecisionIndex;
            size_t mScalingIndex;
            size_t mSilhouetteIndex;
            size_t mEarlyIndex;
            mutable SweepVariant mCached;
        };

        Iterator begin() const
        {
            if (variantCount() == 0)
            {
                return end();
            }

            return Iterator(this, 0, 0, 0, 0);
        }

        Iterator end() const
        {
            return Iterator(this, mPrecisions.size(), 0, 0, 0);
        }
    };

    std::vector<VisibilityExactQueryConfiguration::PrecisionType> allPrecisions()
    {
        return
        {
            VisibilityExactQueryConfiguration::DOUBLE,
            VisibilityExactQueryConfiguration::FLOAT
        #ifdef ENABLE_LEDA
            , VisibilityExactQueryConfiguration::LEDA_REAL
        #endif
        #ifdef ENABLE_GMP
            , VisibilityExactQueryConfiguration::GMP_FLOAT
            , VisibilityExactQueryConfiguration::GMP_RATIONAL
        #endif
        #ifdef ENABLE_MPFR
            , VisibilityExactQueryConfiguration::MPFR
        #endif
        };
    }

    std::vector<VisibilityExactQueryConfiguration::PrecisionType> reducedPrecisions()
    {
        return
        {
            VisibilityExactQueryConfiguration::DOUBLE
        #ifdef ENABLE_MPFR
            , VisibilityExactQueryConfiguration::MPFR
        #endif
        #ifdef ENABLE_GMP
            , VisibilityExactQueryConfiguration::GMP_FLOAT
        #endif
        };
    }

    const std::vector<VisibilityExactQueryConfiguration::SolverType> cSolvers =
    {
        VisibilityExactQueryConfiguration::EXACT_APERTURE_FINDER,
        VisibilityExactQueryConfiguration::EXACT_SEQUENTIAL_SOLVER
    };

    const char* cSolverNames[] = { "aperture", "sequential" };

    std::string listAnswers(const std::vector<VisibilityResult>& answers, const char* lead)
    {
        std::ostringstream line;
        line << lead;

        for (size_t solverIndex = 0; solverIndex < cSolvers.size(); solverIndex++)
        {
            line << (solverIndex == 0 ? " " : ", ") << cSolverNames[solverIndex] << " " << DemoConfiguration::resultName(answers[solverIndex]);
        }

        return line.str();
    }

    std::string validate(const std::vector<VisibilityResult>& answers, VisibilityResult expected,
                         std::vector<size_t>& caseMatches, std::vector<size_t>& precisionMatches)
    {
        bool anyFailure = false;

        for (size_t solverIndex = 0; solverIndex < cSolvers.size(); solverIndex++)
        {
            if (answers[solverIndex] == expected)
            {
                caseMatches[solverIndex]++;
                precisionMatches[solverIndex]++;
            }

            if (answers[solverIndex] == FAILURE)
            {
                anyFailure = true;
            }
        }

        std::ostringstream report;

        if (anyFailure)
        {
            report << "\n    FAIL" << listAnswers(answers, " a solver answered FAILURE:");
        }

        if (expected == UNKNOWN)
        {
            report << "\n    FAIL" << listAnswers(answers, " the case has no expectation,");
        }
        else
        {
            for (size_t solverIndex = 0; solverIndex < cSolvers.size(); solverIndex++)
            {
                if (answers[solverIndex] != expected)
                {
                    report << "\n    FAIL " << cSolverNames[solverIndex] << " answered " << DemoConfiguration::resultName(answers[solverIndex])
                           << ", expected " << DemoConfiguration::resultName(expected);
                }
            }
        }

        return report.str();
    }

    struct CaseOutcome
    {
        VisibilityResult mResult;
        std::string mDescription;
    };

    struct TestResults
    {
        std::vector<CaseOutcome> mOutcomes;
        bool mStopRequested = false;

        void append(VisibilityResult result, const std::string& description)
        {
            CaseOutcome outcome;
            outcome.mResult = result;
            outcome.mDescription = description;

            mOutcomes.push_back(outcome);
        }

        void adopt(TestResults& other)
        {
            for (size_t i = 0; i < other.mOutcomes.size(); i++)
            {
                mOutcomes.push_back(std::move(other.mOutcomes[i]));
            }

            other.mOutcomes.clear();
        }

        size_t failureCount() const
        {
            size_t count = 0;

            for (size_t i = 0; i < mOutcomes.size(); i++)
            {
                if (mOutcomes[i].mResult == FAILURE)
                {
                    count++;
                }
            }

            return count;
        }

        void display() const
        {
            for (size_t i = 0; i < mOutcomes.size(); i++)
            {
                std::cout << std::right << std::setw(cResultWidth) << DemoConfiguration::resultName(mOutcomes[i].mResult) << " | " << mOutcomes[i].mDescription << std::endl;
            }
        }
    };

    void runVariant(const SweepVariant& sweepVariant, const DemoConfiguration& configuration, VisibilityResult expected, const std::string& name, size_t caseIndex,
                    const char* tag, std::vector<size_t>& caseMatches, std::vector<size_t>& precisionMatches, TestResults& results)
    {
        DemoConfiguration settings = configuration;
        settings.globalScaling = sweepVariant.mGlobalScaling;
        settings.silhouetteOptimisation = sweepVariant.mSilhouetteOptimisation;
        settings.detectApertureOnly = sweepVariant.mDetectApertureOnly;
        settings.precisionType = sweepVariant.mPrecision;

        std::ostringstream description;
        description << " gs=" << settings.globalScaling
                    << " silh=" << (settings.silhouetteOptimisation ? "on" : "off")
                    << " early=" << (settings.detectApertureOnly ? "on" : "off")
                    << " prec=" << DemoHelper::precisionTag(settings.precisionType);

        std::vector<float> v0;
        std::vector<float> v1;

        DemoHelper::generatePolygon(v0, settings.vertexCount0, settings.scaling, settings.phi - cSourceHalfTurn, settings.globalScaling);
        DemoHelper::generatePolygon(v1, settings.vertexCount1, settings.scaling, settings.phi, settings.globalScaling);

        HelperTriangleMeshContainer* meshContainer = DemoHelper::createScene(settings.sceneIndex, settings.globalScaling);
        V_ASSERT(meshContainer != nullptr);

        GeometryOccluderSet* occluderSet = DemoHelper::createOccluderSet(meshContainer);

        srand(0);

        std::vector<VisibilityResult> answers(cSolvers.size());
        for (size_t solverIndex = 0; solverIndex < cSolvers.size(); solverIndex++)
        {
            answers[solverIndex] = solveWith(occluderSet, v0, v1, settings, cSolvers[solverIndex]);
        }

        const std::string report = validate(answers, expected, caseMatches, precisionMatches);

        if (report.empty())
        {
            results.append(expected, description.str());
        }
        else
        {
            results.append(FAILURE, description.str() + report);
        }

        delete occluderSet;
        delete meshContainer;
    }

    const int cVariantTimeoutSeconds = 2;

    struct VariantRun
    {
        VariantRun() : mCaseMatches(cSolvers.size(), 0), mPrecisionMatches(cSolvers.size(), 0)
        {
        }

        TestResults mResults;
        std::vector<size_t> mCaseMatches;
        std::vector<size_t> mPrecisionMatches;
    };

    void runVariantWithTimeout(const SweepVariant& sweepVariant, const DemoConfiguration& configuration, VisibilityResult expected, const std::string& name, size_t caseIndex,
                               const char* tag,std::vector<size_t>& caseMatches, std::vector<size_t>& precisionMatches, TestResults& results)
    {
        auto run = std::make_shared<VariantRun>();

        auto finished = std::make_shared<std::promise<void>>();
        std::future<void> signal = finished->get_future();

        std::thread worker([run, finished, sweepVariant, configuration, expected, name, caseIndex, tag]()
        {
            runVariant(sweepVariant, configuration, expected, name, caseIndex, tag,
                       run->mCaseMatches, run->mPrecisionMatches, run->mResults);

            finished->set_value();
        });

        if (signal.wait_for(std::chrono::seconds(cVariantTimeoutSeconds)) == std::future_status::timeout)
        {
            worker.detach();

            std::ostringstream description;
            description << " gs=" << sweepVariant.mGlobalScaling
                        << " silh=" << (sweepVariant.mSilhouetteOptimisation ? "on" : "off")
                        << " early=" << (sweepVariant.mDetectApertureOnly ? "on" : "off")
                        << " prec=" << DemoHelper::precisionTag(sweepVariant.mPrecision)
                        << "\n    FAIL did not answer within " << cVariantTimeoutSeconds << " s, abandoned to the background";

            results.append(FAILURE, description.str());
          
            return;
        }

        worker.join();

        const bool stopRequested = run->mResults.mStopRequested;

        results.adopt(run->mResults);

        if (stopRequested)
        {
            results.mStopRequested = true;
        }

        for (size_t solverIndex = 0; solverIndex < cSolvers.size(); solverIndex++)
        {
            caseMatches[solverIndex] += run->mCaseMatches[solverIndex];
            precisionMatches[solverIndex] += run->mPrecisionMatches[solverIndex];
        }
    }

    void runCase(const std::string& path, const std::string& name, size_t caseIndex, const SweepBounds& sweep, const char* tag, TestResults& results)
    {
        DemoConfiguration configuration;
        configuration.readConfig(path);

        if (configuration.vertexCount0 == 0 || configuration.vertexCount1 == 0 || configuration.sceneIndex >= DemoConfiguration::COUNT)
        {
            results.append(FAILURE, "could not load the geometry of case '" + name + "' from " + path);
            return;
        }

        const VisibilityResult expected = DemoConfiguration::resultFromName(configuration.expected);

        std::vector<size_t> caseMatches(cSolvers.size(), 0);
        std::vector<std::vector<size_t>> precisionMatches(sweep.mPrecisions.size(), std::vector<size_t>(cSolvers.size(), 0));
       
        for (SweepBounds::Iterator it = sweep.begin(); it != sweep.end() && !results.mStopRequested; ++it)
        {
            runVariantWithTimeout(*it, configuration, expected, name, caseIndex,
                                  tag, caseMatches, precisionMatches[it.mPrecisionIndex], results);

        }
    }
}

bool SolverConsistencyTest(std::string& errorMessage)
{
    std::cout << "SolverConsistencyTest" << std::endl;

    const std::string folder = DemoHelper::testCaseFolder();
    const std::vector<std::string> files = DemoHelper::listTestCaseFiles(folder);

    const SweepBounds sweep = { reducedPrecisions(), { 1.0f }, { true, false }, { true, false } };

    TestResults results;

    if (files.empty())
    {
        errorMessage = "no case file found in " + folder;
        std::cout << "    ERROR " << errorMessage << std::endl;
        return false;
    }

    for (size_t caseIndex = 0; caseIndex < files.size(); caseIndex++)
    {
        const std::string& name = files[caseIndex];

        runCase(folder + "/" + name, name, caseIndex, sweep, "solverConsistency", results);
    }

    results.display();

    const size_t failureCount = results.failureCount();

    if (failureCount > 0)
    {
         std::cout << "SolverConsistencyTest FAILED" << std::endl;
        return false;
    }

    std::cout << "SolverConsistencyTest SUCCESS" << std::endl;
    return true;
}

bool SolverEdgeCases(std::string& errorMessage)
{
    std::cout << "SolverEdgeCases" << std::endl;

    const std::string folder = DemoHelper::testCaseFolder() + "/deactivate/on_the_edge";
    const std::vector<std::string> files = DemoHelper::listTestCaseFiles(folder);

    const SweepBounds sweep = { reducedPrecisions(), { 1.0f }, { true, false }, { true, false } };

    TestResults results;

    if (files.empty())
    {
        errorMessage = "no case file found in " + folder;
        std::cout << "    ERROR " << errorMessage << std::endl;
        return false;
    }

    for (size_t caseIndex = 0; caseIndex < files.size(); caseIndex++)
    {
        const std::string& name = files[caseIndex];

        runCase(folder + "/" + name, name, caseIndex, sweep, "solverEdge", results);
    }

    results.display();

    const size_t failureCount = results.failureCount();

    if (failureCount > 0)
    {       
        std::cout << "SolverEdgeCases FAILED" << std::endl;
        return false;
    }

    std::cout << "SolverEdgeCases SUCCESS" << std::endl;
    return true;
}

bool SolverPrecisionTests(std::string& errorMessage)
{
    std::cout << "SolverPrecisionTests" << std::endl;

    const std::string folder = DemoHelper::testCaseFolder() + "/precision_test";
    const std::vector<std::string> files = DemoHelper::listTestCaseFiles(folder);

    const SweepBounds sweep = { reducedPrecisions(), { 1.0f }, { true, false }, { true, false } };

    TestResults results;

    if (files.empty())
    {
        errorMessage = "no case file found in " + folder;
        std::cout << "    ERROR " << errorMessage << std::endl;
        return false;
    }

    for (size_t caseIndex = 0; caseIndex < files.size(); caseIndex++)
    {
        const std::string& name = files[caseIndex];

        runCase(folder + "/" + name, name, caseIndex, sweep, "solverPrecision", results);
    }

    results.display();

    const size_t failureCount = results.failureCount();

    if (failureCount > 0)
    {      
          std::cout << "SolverPrecisionTests FAILED" << std::endl;
        return false;
    }

    std::cout << "SolverPrecisionTests SUCCESS" << std::endl;
    return true;
}
