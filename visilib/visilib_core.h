/*
Visilib, an open source library for exact visibility computation.
Copyright(C) 2020 by Denis Haumont

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

#pragma once

#ifdef EXACT_ARITHMETIC
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/leda_real.h>
typedef CGAL::Lazy_exact_nt<leda_real> exact;
#endif

#include <iostream>
#include <string>
#include <assert.h>

namespace visilib
{
#define V_ASSERT(X) assert(X)

#ifdef LOG_INFORMATION
    inline void V_LOG(const std::string& message, int aDeep = 0)
    {
        for (int i = 0; i < aDeep; i++)
        {
            std::cout << "  ";
        }
        std::cout << message << std::endl;
    }
#else
    inline void V_LOG(std::ostream& out, const std::string& message, std::string prefix = "")
    {
//        out << "["<<prefix<<"] " << message << std::endl;
    }
#endif    
}