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

 if (EXISTS "$ENV{EMBREE_INSTALL_DIR}/include/embree3/rtcore.h" AND EXISTS "$ENV{EMBREE_INSTALL_DIR}/lib/embree3.lib" )
    message("Embree found: "  $ENV{EMBREE_INSTALL_DIR})
    set(EMBREE_INCLUDE_DIRS $ENV{EMBREE_INSTALL_DIR}/include) 
    set(EMBREE_LIBRARY $ENV{EMBREE_INSTALL_DIR}/lib/embree3.lib) 
    set(embree_FOUND True) 
else()
    message("Could not locate embree at " $ENV{EMBREE_INSTALL_DIR})
 
    # under Windows:
    message("Please define the environment variable EMBREE_INSTALL_DIR")
endif()