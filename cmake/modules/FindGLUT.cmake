# Visilib, an open source library for exact visibility computation.
# Copyright(C) 2017 by Denis Haumont
#
# This file is part of Visilib.
#
# Visilib is free software : you can redistribute it and / or modify
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

#ADD_LIBRARY(GLUT_LIBRARY STATIC IMPORTED)
#SET_TARGET_PROPERTIES(GLUT_LIBRARY PROPERTIES
  #                    IMPORTED_LOCATION  ${CMAKE_CURRENT_BINARY_DIR}/../demo/external/freeglut/x64/freeglut.lib)

set(GLUT_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/../demo/external/freeglut/x64/freeglut.lib)
set(GLUT_FOUND true)
