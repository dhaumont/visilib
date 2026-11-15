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

#pragma once

#include <cstddef>

#include "math_vector_3.h"

namespace visilib
{
    /** @brief Compile time description of the cell dimension the Plucker point type lives in.

    The six dimensional space encodes cells of four dimensions: a vertex of a cell lies
    at four facets, and an edge joins two vertices sharing three of them. A line family
    whose members all satisfy the quadric, like MathPluckerParallel6, encodes its cells
    on the two dimensional screen of the family: a vertex of a cell lies at two boundary
    lines, and an edge joins two vertices sharing one. The polytopes read the difference
    here rather than carry the general numbers everywhere.

    cAllLinesReal tells whether every value of the class is a real line, which is what
    makes the quadric corrections of the general space nullities on such a family. */

    template<class P>
    struct MathPluckerCellTraits
    {
        static constexpr size_t cCommonFacetsForEdge = 3;
        static constexpr bool cAllLinesReal = false;
    };

    /** @brief Carry a 3D edge, given by its two endpoints, as a plane of the space of P

    The general point stores the edge itself, its two endpoints being two points of it.
    A line family cannot store a line whose direction is not that of the family as one
    of its rays, and stores the edge as the carrier it imposes on them instead. The
    specialization lives with the family class it knows about. */

    template<class P>
    struct PluckerEdgeCarrier
    {
        static P fromPoints(const MathVector3d& a, const MathVector3d& b)
        {
            return P(a, b);
        }
    };
}
