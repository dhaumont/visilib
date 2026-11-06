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

/* The visitor returns true when a test name should run, which lets a single test be
selected from the command line without rebuilding:

    visilibTest                         run every test
    visilibTest Solver                  run every test whose name contains "Solver"
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
    { "SolverConsistencyTest", &SolverConsistencyTest }
};

static void printUsage(const char* program)
{
    std::cout << "usage: " << program << " [filter] [-l] [-h]" << std::endl
              << "  runs the test suites of the library, one line of result per test" << std::endl
              << "  filter   substring of a test name, e.g. '" << program << " Solver'" << std::endl
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
