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

set(Src
   demo_debug_visualisation_gl.h
   demo_debug_visualisation_gl.cpp
   demo_helper.h
   demo_helper.cpp
   demo_main.cpp
   demo_viewer_glut.h

   external/imgui/imgui_demo.cpp
   external/imgui/imgui.cpp
   external/imgui/misc/cpp/imgui_stdlib.cpp
   external/imgui/imgui_tables.cpp
   external/imgui/backends/imgui_impl_null.cpp
   external/imgui/backends/imgui_impl_opengl2.cpp
   external/imgui/backends/imgui_impl_glut.cpp
   external/imgui/imgui_draw.cpp
   external/imgui/imgui_widgets.cpp
   )

include_directories( external/freeglut/include/ )
include_directories( ../visilib/ )
include_directories( external/imgui )

add_executable( ${VISILIB_DEMO_NAME}  ${Src})
target_link_libraries(${VISILIB_DEMO_NAME} ${LIBS} )

if (USE_MPFR)
  target_link_libraries (${VISILIB_DEMO_NAME} mpfr)
endif()

if (USE_GMP)
  target_link_libraries (${VISILIB_DEMO_NAME} gmp)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_C_FLAGS_DEBUG "-g -O0")

target_compile_options(${VISILIB_DEMO_NAME} PRIVATE -g -rdynamic -fno-omit-frame-pointer)

