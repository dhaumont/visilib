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
    /**@brief The relative position of a point with regards to another geometric entity*/
    enum GeometryPositionType
    {
        ON_BOUNDARY = 1,        /**< @brief The points lies on the geometric entity*/
        ON_NEGATIVE_SIDE = 0,   /**< @brief The points lies on the negative side of geometric entity*/
        ON_POSITIVE_SIDE = 2    /**< @brief The points lies on the postive side of the geometric entity*/
    };
}
