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

#include "visibility_solver.h"
#include "visibility_exact_query.h"

namespace visilib
{
    /** @brief TODO */

    template<class P, class S>
    class VisibilityAggressiveSolver : public VisibilitySolver<P, S>
    {
    public:
        VisibilityAggressiveSolver(VisibilityExactQuery_<P, S>* aSolver,
            S tolerance,
            bool detectApertureOnly);

        VisibilityResult resolve();
    private:
        void extractStabbingLines(PluckerPolyhedron<P>* myPolyhedron, PluckerPolytope<P>* aPolytope);

        bool mDetectApertureOnly;
        S mTolerance;
    };

    template<class P, class S>
    VisibilityAggressiveSolver<P, S>::VisibilityAggressiveSolver(VisibilityExactQuery_<P, S>* mQuery, S tolerance, bool detectApertureOnly)
        : VisibilitySolver<P, S>(mQuery),
        mTolerance(tolerance),
        mDetectApertureOnly(detectApertureOnly)
    {
    }
}
