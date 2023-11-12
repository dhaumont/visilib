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

create_test_sourcelist(Tests
    ./src/test_full.cpp
    ./src/test_math.cpp
    ./src/test_visibility.cpp
    ../demo/demo_helper.cpp
)

include_directories( ../visilib/ )

add_executable(visilibTest ${Tests})

target_link_libraries(visilibTest ${LIBS})

#add_test(NAME stc_test_math COMMAND ${KIT}CxxTests stc_test_math})
add_test(NAME visilibTest COMMAND visilibTest)

# remove the test driver source file
set (TestsToRun ${Tests})
list ( REMOVE_ITEM TestsToRun src/test_full.cxx)

# Add all the ADD_TEST for each test
foreach (test ${TestsToRun})
    get_filename_component (TName ${test} NAME_WE)
    add_test (NAME ${TName} COMMAND test_full ${TName})
endforeach ()
