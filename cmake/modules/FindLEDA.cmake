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

if (EXISTS "$ENV{LEDA_DIR}/incl/LEDA/numbers/real.h" AND EXISTS "$ENV{LEDA_DIR}/leda_mt.lib" )
    message("LEDA found: "  $ENV{LEDA_DIR})
 
    set (LEDA_INCLUDE_DIR $ENV{LEDA_DIR}/incl)
    set (LEDA_LIBRARY $ENV{LEDA_DIR}/leda_mt.lib)
    set(LEDA_FOUND True)
else()
    message("Could not locate LEDA at " $ENV{LEDA_DIR} )
    message("Please define the environment variable LEDA_DIR")
endif()

