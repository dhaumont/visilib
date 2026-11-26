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

namespace visilib
{
    /** @brief How the Plucker points are brought to a common scale (see MathPluckerNormalization).

    A Plucker point is defined up to a positive factor: the choice changes no predicate, only the magnitudes the
    tolerance is compared to and the rounding of the constructions.*/
    enum NormalizationType
    {
        NO_NORMALIZATION = 0,           /**< @brief The points keep the scale they are built with*/
        HYPERSPHERE_NORMALIZATION = 1,  /**< @brief The points are projected on the unit hypersphere (a square root)*/
        CHART_NORMALIZATION = 2         /**< @brief The cell points are divided by a linear form positive on every line from the source to the receiver, the hyperplanes by a power of two (no square root, rational)*/
    };
}
