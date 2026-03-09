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

#include "geometry_convex_polygon.h"
#include "silhouette_mesh_face.h"
#include "math_plucker_2.h"
#include "math_predicates.h"
#include "plucker_polyhedron.h"
#include "plucker_polytope.h"
#include "plucker_polytope_complex.h"
#include "plucker_polytope_splitter.h"
#include "visibility_exact_query.h"
#include "silhouette.h"
#include <fstream>
#include <sstream>
#include <map>
namespace visilib
{
    template<class P>
    class  PluckerPolytope;

    /** @brief  */

    template<class P, class S>
    class VisibilitySequentialSolver : public VisibilitySolver<P, S>
    {
    public:
        VisibilitySequentialSolver(VisibilityExactQuery_<P, S>* aSolver,
            bool normalization,
            S tolerance,
            bool detectApertureOnly);

        VisibilityResult resolve();
    private:
      
        bool mNormalization;
        bool mDetectApertureOnly;
        S mTolerance;
    };

    template<class P, class S>
    VisibilitySequentialSolver<P, S>::VisibilitySequentialSolver(VisibilityExactQuery_<P, S>* mQuery, bool normalization, S tolerance, bool detectApertureOnly)
        : VisibilitySolver<P, S>(mQuery),
        mNormalization(normalization),
        mTolerance(tolerance),
        mDetectApertureOnly(detectApertureOnly)
    {
    }

    template <class P, class S>
    VisibilityResult VisibilitySequentialSolver<P, S>::resolve()
    {
        VisibilityResult myGlobalResult = UNKNOWN;

       return myGlobalResult;
    }
}