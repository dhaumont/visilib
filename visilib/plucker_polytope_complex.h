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

#include <stack>
#include "math_plucker_6.h"

namespace visilib
{
    template<class S>
    class PluckerPolytope;
    template<class S>
    class PluckerPolyhedron;
    class HelperStatisticCollector;
    template<class P, class S>
    class VisibilityExactQuery_;
    class GeometryConvexPolygon;

   
    /** @brief Represents a complex of polytopes in Plucker space.

    Currently only the initial polytope and the polyhedron in Pluker space are stored explicitely
    */

    template<class P>
    class PluckerPolytopeComplex
    {
    public:
        PluckerPolytopeComplex();

        ~PluckerPolytopeComplex();

        bool hasPolytope()const
        {
            return mRoot != nullptr;
        }

        void setRoot(PluckerPolytope<P>* polytope)
        {
            V_ASSERT(!hasPolytope());
            mRoot = polytope;
        }

        PluckerPolyhedron<P>* getPolyhedron() { return mPolyhedron; }
        PluckerPolytope<P>* getRoot() { return mRoot; }
    private:

        PluckerPolyhedron<P>* mPolyhedron;
        PluckerPolytope<P>* mRoot;
    };

    template<class P>
    inline PluckerPolytopeComplex<P>::PluckerPolytopeComplex()
    {
        mRoot = nullptr;
        mPolyhedron = new PluckerPolyhedron<P>();
    }

    template<class P>
    inline PluckerPolytopeComplex<P>::~PluckerPolytopeComplex()
    {
        delete mRoot;
        delete mPolyhedron;
    }
}
