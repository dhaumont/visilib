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

#include "test_definition.h"

/* The first argument selects the suites by name, the second narrows the sweeps of the solver
suites to the variant numbers their output prints:

    visilibTest                         run every test
    visilibTest Solver                  run every test whose name contains "Solver"
    visilibTest Solver 12,14,15         of those, run only the variants 12, 14 and 15
    visilibTest -l                      list the available tests
    visilibTest -h                      print this help
*/

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
    std::cout << "usage: " << program << " [filter] [ids] [-l] [-h]" << std::endl
              << "  runs the test suites of the library, one line of result per test" << std::endl
              << "  filter   substring of a test name, e.g. '" << program << " Solver'" << std::endl
              << "  ids      comma separated variant numbers to run, e.g. '12,14,15'," << std::endl
              << "           as printed in the '#n' column of the solver suites; every" << std::endl
              << "           variant runs when it is left out" << std::endl
              << "  -l       list the available tests and exit" << std::endl
              << "  -h       print this message and exit" << std::endl << "  tests:" << std::endl;

    for (size_t i = 0; i < sizeof(cTests) / sizeof(cTests[0]); i++)
    {
        std::cout << "    " << cTests[i].mName << std::endl;
    }
}

int main(int argc, char** argv)
{
    if (argc > 1 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help"))
    {
        printUsage(argv[0]);
        return 0;
    }

    if (argc > 1 && std::string(argv[1]) == "-l")
    {
        for (size_t i = 0; i < sizeof(cTests) / sizeof(cTests[0]); i++)
        {
            std::cout << cTests[i].mName << std::endl;
        }
        return 0;
    }

    const std::string filter = argc > 1 ? argv[1] : std::string();

    /* A second argument narrows the solver sweeps to the variant numbers they print. */
    if (argc > 2)
    {
        selectSolverTestVariants(argv[2]);
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
