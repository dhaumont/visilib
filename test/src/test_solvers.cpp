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
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <iterator>
#include <limits>
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
    const size_t cIdWidth = 6;
    const size_t cResultWidth = 8; 
    const size_t cAnswerWidth = 18;   // answer plus the optional " != expected" part
    const float cSourceHalfTurn = 3.14519f;

    /* The variants of the solver sweeps are numbered, from 1, in the order they would run, and
    the command line can name a few of those numbers to run only them. A variant that is not
    named still spends its number, which is what keeps a number meaning the same variant from one
    run to the next; naming no number runs everything, which is the default. */
    struct VariantSelection
    {
        std::vector<size_t> mRequested;   // empty means every variant
        std::vector<bool> mSeen;          // mSeen[id], ids start at 1
        size_t mHandedOut = 0;            // how many ids the sweeps have spent
        size_t mExecuted = 0;             // how many of those have run
    };

    VariantSelection gSelection;

    size_t nextVariantId()
    {
        const size_t id = ++gSelection.mHandedOut;

        if (gSelection.mSeen.size() <= id)
        {
            gSelection.mSeen.resize(id + 1, false);
        }

        return id;
    }

    bool isVariantSelected(size_t id)
    {
        if (gSelection.mRequested.empty())
        {
            return true;
        }

        const bool selected = std::find(gSelection.mRequested.begin(), gSelection.mRequested.end(), id) != gSelection.mRequested.end();

        if (selected)
        {
            gSelection.mSeen[id] = true;
        }

        return selected;
    }

    /* What a selection leaves out is worth stating, because the run looks short otherwise, and
    a number that no variant carried says the sweep has changed since the number was read. */
    void displaySelectionReport()
    {
        if (gSelection.mRequested.empty())
        {
            return;
        }

        std::cout << "    " << gSelection.mExecuted << " of " << gSelection.mHandedOut << " variants selected" << std::endl;

        std::ostringstream missing;

        for (size_t i = 0; i < gSelection.mRequested.size(); i++)
        {
            const size_t id = gSelection.mRequested[i];

            if (id >= gSelection.mSeen.size() || !gSelection.mSeen[id])
            {
                missing << (missing.str().empty() ? "" : ", ") << id;
            }
        }

        if (!missing.str().empty())
        {
            std::cout << "    no variant with id " << missing.str() << " was run" << std::endl;
        }
    }
    
    VisibilityResult solveWith(GeometryOccluderSet* occluderSet, const std::vector<float>& v0, const std::vector<float>& v1,
                               const DemoConfiguration& configuration, VisibilityExactQueryConfiguration::SolverType solverType,
                               double timeoutSeconds, bool& budgetExhausted)
    {
        VisibilityExactQueryConfiguration query;
        query.detectApertureOnly = configuration.detectApertureOnly;
        query.hyperSphereNormalization = configuration.normalization;
        query.silhouetteOptimization = configuration.silhouetteOptimisation;
        query.precision = configuration.precisionType;
        query.tolerance = configuration.tolerance;
        query.confidenceValue = configuration.confidenceValue;
        query.solverType = solverType;
        query.displaySplitOccurence = std::numeric_limits<int>::max();

        const double side = 2.0 * configuration.scaling;
        const double aperture = configuration.minimumApertureSize;
        query.minimumNormalizedApertureSize = (aperture * aperture) / (side * side);
 
        HelperVisualDebugger debugger;
        debugger.getTimeout().setMaximumDuration(timeoutSeconds);
        debugger.setDisplayStatistics(false);

        const VisibilityResult result = visilib::areVisible(occluderSet, &v0[0], v0.size() / 3, &v1[0], v1.size() / 3, query, &debugger);

         budgetExhausted = debugger.getTimeout().isOver();

        return result;
    }

    struct SweepVariant
    {
        VisibilityExactQueryConfiguration::PrecisionType mPrecision = VisibilityExactQueryConfiguration::DOUBLE;
        float mGlobalScaling = 1.0f;
        bool mSilhouetteOptimisation = false;
        bool mDetectApertureOnly = false;
        VisibilityExactQueryConfiguration::SolverType mSolver = VisibilityExactQueryConfiguration::EXACT_APERTURE_FINDER;
    };

    struct SweepBounds
    {
        std::vector<VisibilityExactQueryConfiguration::PrecisionType> mPrecisions;
        std::vector<float> mGlobalScalings;
        std::vector<bool> mSilhouetteOptimisations;
        std::vector<bool> mDetectApertureOnly;
        std::vector<VisibilityExactQueryConfiguration::SolverType> mSolvers;

        size_t variantCount() const
        {
            return mGlobalScalings.size() * mSilhouetteOptimisations.size() * mDetectApertureOnly.size() * mSolvers.size()* mPrecisions.size();
        }
  
        struct Iterator
        {
            typedef std::forward_iterator_tag iterator_category;
            typedef SweepVariant value_type;
            typedef std::ptrdiff_t difference_type;
            typedef const SweepVariant* pointer;
            typedef const SweepVariant& reference;

            Iterator() : mSweep(nullptr), mPrecisionIndex(0), mScalingIndex(0), mSilhouetteIndex(0), mEarlyIndex(0), mSolverIndex(0)
            {
            }

            Iterator(const SweepBounds* sweep, size_t precisionIndex, size_t scalingIndex, size_t silhouetteIndex, size_t earlyIndex, size_t solverIndex)
                : mSweep(sweep), mPrecisionIndex(precisionIndex), mScalingIndex(scalingIndex), mSilhouetteIndex(silhouetteIndex), mEarlyIndex(earlyIndex), mSolverIndex(solverIndex)
            {
            }

            SweepVariant at(size_t precisionIndex, size_t scalingIndex, size_t silhouetteIndex, size_t earlyIndex, size_t solverIndex) const
            {
                SweepVariant variant;
                variant.mPrecision = mSweep->mPrecisions[precisionIndex];
                variant.mGlobalScaling = mSweep->mGlobalScalings[scalingIndex];
                variant.mSilhouetteOptimisation = mSweep->mSilhouetteOptimisations[silhouetteIndex];
                variant.mDetectApertureOnly = mSweep->mDetectApertureOnly[earlyIndex];
                variant.mSolver = mSweep->mSolvers[solverIndex];

                return variant;
            }

            reference operator*() const
            {
                mCached = at(mPrecisionIndex, mScalingIndex, mSilhouetteIndex, mEarlyIndex, mSolverIndex);
                return mCached;
            }

            pointer operator->() const
            {
                mCached = at(mPrecisionIndex, mScalingIndex, mSilhouetteIndex, mEarlyIndex, mSolverIndex);
                return &mCached;
            }

            Iterator& operator++()
            {
                mSolverIndex++;

                if (mSolverIndex == mSweep->mSolvers.size())
                {
                    mSolverIndex = 0;
                    mEarlyIndex++;
                }

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
                       mSilhouetteIndex == other.mSilhouetteIndex && mEarlyIndex == other.mEarlyIndex && mSolverIndex == other.mSolverIndex;
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
            size_t mSolverIndex;
            mutable SweepVariant mCached;
        };

        Iterator begin() const
        {
            if (variantCount() == 0)
            {
                return end();
            }

            return Iterator(this, 0, 0, 0, 0, 0);
        }

        Iterator end() const
        {
            return Iterator(this, mPrecisions.size(), 0, 0, 0, 0);
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

    const char* solverName(VisibilityExactQueryConfiguration::SolverType solver)
    {
        switch (solver)
        {
            case VisibilityExactQueryConfiguration::EXACT_SEQUENTIAL_SOLVER: return "sequential";
            default: return "aperture";
        }
    }

    std::string variantDescription(const SweepVariant& sweepVariant)
    {
        std::ostringstream description;
        description << " [ scaling: " << sweepVariant.mGlobalScaling
                    << ", silhouette: " << (sweepVariant.mSilhouetteOptimisation ? "ON" : "OFF")
                    << ", earlyStop: " << (sweepVariant.mDetectApertureOnly ? "ON" : "OFF")
                    << ", precision: " << DemoHelper::precisionTag(sweepVariant.mPrecision)
                    << ", solver: " << solverName(sweepVariant.mSolver)
                    << " ]";

        return description.str();
    }
    enum TestResult
    {
        TEST_RESULT_PASS,
        TEST_RESULT_FAIL,
        TEST_RESULT_UNKNOWN,
        TEST_RESULT_TIMEOUT
    };

    struct CaseOutcome
    {
        size_t mId;
        TestResult  mResult;
        VisibilityResult mAnswer;
        VisibilityResult mExpected;
        std::string mDescription;

        CaseOutcome(size_t id, TestResult result, VisibilityResult answer, VisibilityResult expected, const std::string& description)
         : mId(id), mResult(result), mAnswer(answer), mExpected(expected), mDescription(description)
        {
        }

        const std::string toStr(TestResult p) const
        {
            switch (p) 
            {
            case TEST_RESULT_PASS:
                return "PASS";
                break;
            case TEST_RESULT_FAIL:
                return "FAIL";
                break;
            case TEST_RESULT_UNKNOWN:
                return "UNKNOWN";
                break;
            case TEST_RESULT_TIMEOUT:
                return "TIMEOUT";
                break;
            }
            return "UNKNOWN";
        }
        /* The line is built once so that the pending and the final display share the
        exact same fixed width layout, which lets display() overwrite the pending
        line cleanly with a carriage return. */
        std::string line() const
        {
            std::ostringstream answer;
            answer << DemoConfiguration::toStr(mAnswer);

            if (mResult == TEST_RESULT_FAIL)
            {
                answer << " != " << DemoConfiguration::toStr(mExpected);
            }

            const std::string id = "#" + std::to_string(mId);

            std::ostringstream out;
            out   << " " << std::right << std::setw(cIdWidth) << id << " | " << std::left
                  << std::setw(cResultWidth) << toStr(mResult) << " | "
                  << std::setw(cAnswerWidth) << answer.str() << " | "
                  << mDescription;

            return out.str();
        }

        void display() const
        {
            // The leading carriage return replaces the pending line when one was left
            // unterminated on the current line, and is a no-op on a fresh line.
            std::cout << "\r" << line() << std::endl;
        }

        /* Outputs the line without terminating it, so that the subsequent display()
        call overwrites it once the answer is known. */
        void displayPending() const
        {
            std::cout << line() << std::flush;
        }
    };
    
    struct TestResults
    {
        std::vector<CaseOutcome> mOutcomes;
        bool mStopRequested = false;
        
        void append(const CaseOutcome& outcome)
        {
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
                if (mOutcomes[i].mResult == TEST_RESULT_FAIL)
                {
                    count++;
                }
            }

            return count;
        }

 
        void display() const
        {
            for (const auto& outcome : mOutcomes)
            {
                outcome.display();
            }
        }
    };

    void runVariant(size_t variantId, const SweepVariant& sweepVariant, const DemoConfiguration& configuration, VisibilityResult expected, TestResults& results)
    {
        DemoConfiguration settings = configuration;
        settings.globalScaling = sweepVariant.mGlobalScaling;
        settings.silhouetteOptimisation = sweepVariant.mSilhouetteOptimisation;
        settings.detectApertureOnly = sweepVariant.mDetectApertureOnly;
        settings.precisionType = sweepVariant.mPrecision;

        const std::string description = variantDescription(sweepVariant);

        std::vector<float> v0;
        std::vector<float> v1;

        DemoHelper::generatePolygon(v0, settings.vertexCount0, settings.scaling, settings.phi - cSourceHalfTurn, settings.globalScaling);
        DemoHelper::generatePolygon(v1, settings.vertexCount1, settings.scaling, settings.phi, settings.globalScaling);

        HelperTriangleMeshContainer* meshContainer = DemoHelper::createScene(settings.sceneIndex, settings.globalScaling);
        V_ASSERT(meshContainer != nullptr);

        GeometryOccluderSet* occluderSet = DemoHelper::createOccluderSet(meshContainer);

        srand(0);

        CaseOutcome outcome(variantId, TEST_RESULT_UNKNOWN, UNKNOWN, expected, description);
        
        outcome.displayPending();

         bool budgetExhausted = false;

        const auto start = std::chrono::steady_clock::now();
        outcome.mAnswer = solveWith(occluderSet, v0, v1, settings, sweepVariant.mSolver, queryBudgetSeconds(), budgetExhausted);
        const long elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();

        outcome.mDescription += " (" + std::to_string(elapsedMs) + "ms)";

        if (expected != UNKNOWN)
        {
            if (outcome.mAnswer == FAILURE && budgetExhausted)
            {
               outcome.mResult = TEST_RESULT_TIMEOUT;
            }
            else
            {
                outcome.mResult = (outcome.mAnswer == expected) ? TEST_RESULT_PASS : TEST_RESULT_FAIL;
            }
        }
     
        outcome.display();   
        results.append(outcome);
  
        delete occluderSet;
        delete meshContainer;
    }

    void runCase(const std::string& path, const std::string& name, const SweepBounds& sweep, TestResults& results)
    {
        DemoConfiguration configuration;
        configuration.readConfig(path);

        if (configuration.vertexCount0 == 0 || configuration.vertexCount1 == 0 || configuration.sceneIndex >= DemoConfiguration::COUNT)
        {
            std::cout << "Case '" << name <<  std::endl;

            const CaseOutcome outcome(nextVariantId(), TEST_RESULT_FAIL, UNKNOWN, UNKNOWN, "Could not load the geometry of case '" + name + "' from " + path);
            results.append(outcome);
            outcome.display();
            return;
        }

        const VisibilityResult expected = DemoConfiguration::resultFromName(configuration.expected);

        bool headerDisplayed = false;

        for (SweepBounds::Iterator it = sweep.begin(); it != sweep.end() && !results.mStopRequested; ++it)
        {
            const size_t variantId = nextVariantId();

            if (!isVariantSelected(variantId))
            {
                continue;
            }

            if (!headerDisplayed)
            {
                std::cout << "Case '" << name <<  std::endl;
                headerDisplayed = true;
            }

            gSelection.mExecuted++;

            runVariant(variantId, *it, configuration, expected, results);
        }
    }
}

void selectSolverTestVariants(const std::string& idList)
{
    std::istringstream stream(idList);
    std::string token;

    while (std::getline(stream, token, ','))
    {
        token.erase(std::remove_if(token.begin(), token.end(), [](char c) { return c == ' ' || c == '\t'; }), token.end());

        if (token.empty())
        {
            continue;
        }

        std::istringstream number(token);
        size_t id = 0;

        if (!(number >> id) || id == 0)
        {
            std::cerr << "    ignoring '" << token << "': not a variant number" << std::endl;
            continue;
        }

        gSelection.mRequested.push_back(id);
    }
}

/* The budget lives here rather than in the header constant it starts from, because the runner
changes it from the command line and both suites read it back. The unlimited value is the one
HelperTimeout already carries when nobody sets a budget on it. */
double gQueryBudgetSeconds = cDefaultQueryBudgetSeconds;

double queryBudgetSeconds()
{
    return gQueryBudgetSeconds;
}

void disableQueryBudget()
{
    gQueryBudgetSeconds = std::numeric_limits<double>::max();
}

bool SolverConsistencyTest(std::string& errorMessage)
{
    std::cout << "SolverConsistencyTest" << std::endl;

    const std::string folder = DemoHelper::testCaseFolder();
    const std::vector<std::string> files = DemoHelper::listTestCaseFiles(folder);

    const SweepBounds sweep = { reducedPrecisions(), { 1.0f }, { true, false }, { true, false }, cSolvers };

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

        runCase(folder + "/" + name, name, sweep, results);
    }

   
    displaySelectionReport();

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

    const SweepBounds sweep = { reducedPrecisions(), { 1.0f }, { true, false }, { true, false }, cSolvers };

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

        runCase(folder + "/" + name, name, sweep, results);
    }

    results.display();

    displaySelectionReport();

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

    const SweepBounds sweep = { reducedPrecisions(), { 1.0f }, { true, false }, { true, false }, cSolvers };

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

        runCase(folder + "/" + name, name, sweep, results);
    }

    results.display();

    displaySelectionReport();

    const size_t failureCount = results.failureCount();

    if (failureCount > 0)
    {      
          std::cout << "SolverPrecisionTests FAILED" << std::endl;
        return false;
    }

    std::cout << "SolverPrecisionTests SUCCESS" << std::endl;
    return true;
}
