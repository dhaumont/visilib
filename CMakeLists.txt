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


# top-level CMakeLists.txt
# !! To be run from 'cmakeout' subdirectory !!

cmake_minimum_required (VERSION 3.1 FATAL_ERROR)

set( VISILIB_LIB_NAME "visilib" )
set( VISILIB_DEMO_NAME "visilibDemo" )

set( VERSION_MAJOR 1 )
set( VERSION_MINOR 0 )
set( VERSION_REVISION 0 )

project( ${VISILIB_LIB_NAME} )
project( ${VISILIB_DEMO_NAME} )

if(WIN32)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
endif()


#Set USE_EMBREE flag to True to enable fast raytracing using Intel Embree library (default: false)
set (USE_EMBREE False)

#Set USE_EXACT_ARITHMETIC flag to True to enable exact arithmetic computation with CGAL library (default:false) 
set (USE_EXACT_ARITHMETIC False)

#Comment the following line to disable the sanity integrity checks performed during computations in debug mode (default: checks are activated)
add_definitions(-DINTEGRITY_CHECK)

#Uncomment the following line to perform full validity checks during computations in debug mode (default: checks are not activated)
#add_definitions(-DFULL_VALIDITY_CHECK)

#uncomment the following line to output log information into a debug file during computation (debug and release mode) (default: output is not activated)
#add_definitions(-DOUTPUT_DEBUG_FILE)

#GLUT library (used for the demo part only)
find_package(GLUT)
if (GLUT_FOUND)
    add_definitions(-DUSE_GLUT)
    include_directories(${GLUT_INCLUDE_DIRS})
    link_directories(${GLUT_LIBRARY_DIRS})
    add_definitions(${GLUT_DEFINITIONS})
    set (LIBS ${GLUT_LIBRARIES})
else()
    message(ERROR " GLUT not found!")
endif()

#optional usage of CGAL library
if (USE_EXACT_ARITHMETIC)
    #CGAL library
    find_package(CGAL)
    if (CGAL_FOUND)
        include_directories ( ${CGAL_INCLUDE_DIRS})
        add_definitions(-DCGAL_HEADER_ONLY -DCGAL_NO_GMP)
    endif()

    #BOOST library
    set(Boost_USE_STATIC_LIBS        OFF) # only find static libs
    find_package(Boost)

    if (Boost_FOUND)
        include_directories(${Boost_INCLUDE_DIRS})
    endif()

    find_package(LEDA)
    if (LEDA_FOUND)
        INCLUDE_DIRECTORIES(${LEDA_INCLUDE_DIR})
        set(LIBS ${LIBS} ${LEDA_LIBRARY})
    endif()

    if (Boost_FOUND AND LEDA_FOUND AND CGAL_FOUND)
        add_definitions(-DEXACT_ARITHMETIC)
        add_definitions(-DLEDA_DLL)
        add_definitions(-DCGAL_USE_LEDA)
    endif()
endif()

# optional usage of EMBREE library
if (USE_EMBREE)
    find_package(embree 3.0 REQUIRED)
    if (embree_FOUND)
        INCLUDE_DIRECTORIES(${EMBREE_INCLUDE_DIRS})
        set(LIBS ${LIBS} ${EMBREE_LIBRARY})
        add_definitions(-DEMBREE)
   endif()
endif()

if(NOT WIN32)
    find_package(OpenGL COMPONENTS OpenGL)
    include_directories(${OpenGL_INCLUDE_DIRS})
    link_directories(${OpenGL_LIBRARY_DIRS})
    add_definitions(${OpenGL_DEFINITIONS})
   
   
    set(LIBS ${LIBS} ${OPENGL_LIBRARIES})
    if(NOT OPENGL_FOUND)
        message(ERROR " OPENGL not found!")
    endif(NOT OPENGL_FOUND)

endif(NOT WIN32)

ENABLE_TESTING()


add_subdirectory(visilib)
add_subdirectory(demo)
add_subdirectory(test)

