# Visilib, an open source library for exact visibility computation.
# Copyright(C) 2021 by Denis Haumont
#
# This file is part of Visilib.
#
# Visilib is free software : you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Visilib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Visilib. If not, see <http://www.gnu.org/licenses/>

set(KIT visilibTest)

# ${KIT}CxxTests.cpp

set(Tests
    ./src/test_full.cpp
    ./src/test_math.cpp
    ./src/test_visibility.cpp
    ./src/test_solvers.cpp
    ../demo/demo_helper.cpp
)

include_directories( ../visilib/ )

add_definitions(-DTEST_CASES_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_cases")

add_executable(visilibTest ${Tests})

target_link_libraries(visilibTest ${LIBS})

if (USE_MPFR)
  target_link_libraries (visilibTest mpfr)
endif()

if (USE_GMP)
  target_link_libraries (visilibTest gmp)
endif()
#add_test(NAME stc_test_math COMMAND ${KIT}CxxTests stc_test_math})
add_test(NAME visilibTest COMMAND visilibTest)

# visilibTest is a single binary running the whole suite (see src/test_full.cpp):
# there is no per file executable, so no per file ADD_TEST is registered here.
# Run a subset by restricting the case tables in the test sources, and use
# CTest's -R visilibTest to select the suite.
