# Visilib, an open source library for exact visibility computation.
# Copyright(C) 2020 by Denis Haumont
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

create_test_sourcelist(Tests ${KIT}CxxTests.cpp
    src/test_math.cpp
    )

include_directories( ../visilib/ )

add_executable(visilibTest src/test_math.cpp)
add_executable(visilibTest src/test_visibility.cpp)
add_executable(visilibTest src/test_full.cpp)

add_executable(visilibTest ../demo/demo_helper.cpp)

target_link_libraries(visilibTest ${LIBS})

#add_test(NAME stc_test_math COMMAND ${KIT}CxxTests stc_test_math})
add_test(NAME visilibTest COMMAND visilibTest)
