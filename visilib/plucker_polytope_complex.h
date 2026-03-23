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


        void appendPolytope(PluckerPolytope<P>* polytope)
        {            
            if (polytope == NULL)
            {
                mNeedRepack = true;
            }
            mPolytopes.push_back(polytope);
        }

        PluckerPolyhedron<P>* getPolyhedron() { return mPolyhedron; }
        PluckerPolytope<P>* getRoot() { return mPolytopes[0]; }
        size_t getPolytopeCount() { return mPolytopes.size();}
        PluckerPolytope<P>* getPolytope(size_t index, bool test_nullptr = true) 
        {
            V_ASSERT(!test_nullptr || mPolytopes[index] != NULL);
             return mPolytopes[index];
        }
        void removeLast()
         {
            if (mPolytopes.empty())
            {
                return;
            }
            PluckerPolytope<P>* lastPolytope =  mPolytopes[mPolytopes.size()-1];
            if (lastPolytope != NULL)
            {
                delete lastPolytope;
            }
             mPolytopes.resize(mPolytopes.size()-1);
        }
        void setPolytope(size_t index, PluckerPolytope<P>* p) 
        {
            if (p == NULL)
            {
                mNeedRepack = true;
            }
            mPolytopes[index] = p;
        }
        void repack()
        {
            if (mNeedRepack)
            {
                mNeedRepack = false;
                mPolytopes.erase(std::remove_if(mPolytopes.begin(), mPolytopes.end(), [](PluckerPolytope<P>*  p) { return p == NULL; }), mPolytopes.end());
            }
        }
    private:

        PluckerPolyhedron<P>* mPolyhedron;
        std::vector<PluckerPolytope<P>*> mPolytopes;
        bool mNeedRepack;
    };

    template<class P>
    inline PluckerPolytopeComplex<P>::PluckerPolytopeComplex()
    {
        mPolyhedron = new PluckerPolyhedron<P>();
        mNeedRepack = false;
    }

    template<class P>
    inline PluckerPolytopeComplex<P>::~PluckerPolytopeComplex()
    {
        for (auto myPolytope:mPolytopes)
        {
          delete myPolytope;
        }
        delete mPolyhedron;
    }
}

