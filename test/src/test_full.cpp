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
#include <vector>
#include <iostream>

#include "test_definition.h"

/* The arguments select the suites by name, then the variant numbers of the solver sweeps, and
NO_TIMEOUT lifts the time budget:

    visilibTest                         run every test
    visilibTest Solver                  run every test whose name contains "Solver"
    visilibTest Solver 12,14,15         of those, run only the variants 12, 14 and 15
    visilibTest Solver 12,14 NO_TIMEOUT of those, run them with no time budget
    visilibTest -l                      list the available tests
    visilibTest -h                      print this help
*/

/* Whether a run gives the queries a time budget is not a choice about the tests, so it is said
with its own argument rather than by putting a number where the suite name goes:

    visilibTest SolverConsistencyTest NO_TIMEOUT     answer every variant, however long it takes

The budget, 2.0s, stops a variant that has stopped being interesting; NO_TIMEOUT lifts it, which
is how a variant reported as TIMEOUT gets looked at. */
const char* const cNoTimeoutArgument = "NO_TIMEOUT";

typedef bool (*TestFunction)(std::string&);

struct TestDefinition
{
    const char* mName;
    TestFunction mFunction;
};

const TestDefinition cTests[] =
{
    { "MathCombinatorialTest", &MathCombinatorialTest },
    { "VisibilityTest", &VisibilityTest },
    { "SolverConsistencyTest", &SolverConsistencyTest },
    { "SolverEdgeCases", &SolverEdgeCases },
    { "SolverPrecisionTests", &SolverPrecisionTests }
};

static void printUsage(const char* program)
{
    std::cout << "usage: " << program << " [filter] [ids] [NO_TIMEOUT] [-l] [-h]" << std::endl
              << "  runs the test suites of the library, one line of result per test" << std::endl
              << "  filter   substring of a test name, e.g. '" << program << " Solver'" << std::endl
              << "  ids      comma separated variant numbers to run, e.g. '12,14,15'," << std::endl
              << "           as printed in the '#n' column of the solver suites; every" << std::endl
              << "           variant runs when it is left out" << std::endl
              << "  NO_TIMEOUT" << std::endl
              << "           run the queries without a time budget, so that a variant which" << std::endl
              << "           reports TIMEOUT is answered instead of cut short" << std::endl
              << "  -l       list the available tests and exit" << std::endl
              << "  -h       print this message and exit" << std::endl << "  tests:" << std::endl;

    for (size_t i = 0; i < sizeof(cTests) / sizeof(cTests[0]); i++)
    {
        std::cout << "    " << cTests[i].mName << std::endl;
    }
}

int main(int argc, char** argv)
{
    /* The arguments are the suite name filter and then the variant numbers, in either order with
    NO_TIMEOUT wherever it suits. */
    std::vector<std::string> arguments;
    bool noTimeout = false;

    for (int i = 1; i < argc; i++)
    {
        if (std::string(argv[i]) == cNoTimeoutArgument)
        {
            disableQueryBudget();
            noTimeout = true;
            continue;
        }

        arguments.push_back(argv[i]);
    }

    /* Said once, so a run that takes its time is not read as a run that hangs. */
    if (noTimeout)
    {
        std::cout << "NO_TIMEOUT: the queries run without a time budget" << std::endl;
    }

    const bool hasArgument = !arguments.empty();

    if (hasArgument && (arguments[0] == "-h" || arguments[0] == "--help"))
    {
        printUsage(argv[0]);
        return 0;
    }

    if (hasArgument && arguments[0] == "-l")
    {
        for (size_t i = 0; i < sizeof(cTests) / sizeof(cTests[0]); i++)
        {
            std::cout << cTests[i].mName << std::endl;
        }
        return 0;
    }

    const std::string filter = hasArgument ? arguments[0] : std::string();

    if (arguments.size() > 1)
    {
        selectSolverTestVariants(arguments[1]);
    }

    size_t selectedCount = 0;
    size_t failureCount = 0;

    for (size_t i = 0; i < sizeof(cTests) / sizeof(cTests[0]); i++)
    {
        if (!filter.empty() && std::string(cTests[i].mName).find(filter) == std::string::npos)
        {
            continue;
        }

        selectedCount++;

        std::string errorMessage;

        if (!cTests[i].mFunction(errorMessage))
        {
            std::cout << cTests[i].mName << " FAILED";

            if (!errorMessage.empty())
            {
                std::cout << ": " << errorMessage;
            }

            std::cout << std::endl;
            failureCount++;
        }
    }

    if (selectedCount == 0)
    {
        std::cout << "no test matches '" << filter << "'" << std::endl;
        printUsage(argv[0]);
        return 1;
    }

    std::cout << failureCount << " of " << selectedCount << " tests failed" << std::endl;
    return failureCount > 0 ? 1 : 0;
}
