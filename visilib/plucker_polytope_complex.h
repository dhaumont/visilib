/*
Visilib, an open source library for exact visibility computation.
Copyright(C) 2017 by Denis Haumont

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
    class VisibilityExactQuery;
    class GeometryConvexPolygon;

    /** @brief Represents one node of a BSP tree, the complete tree encodes an occlusion tree in Plucker Space 
    
    Each node stores a polytope, a splitting plane and a pointer to a left and a right polytope
    */
    template<class P>
    class PluckerTreeNode
    {
    public:
        PluckerTreeNode(PluckerPolytope<P>* aPolytope, int aSplittingPlane)
            : mPolytope(aPolytope),
            mSplittingPlane(aSplittingPlane),
            mLeft(nullptr),
            mRight(nullptr),
            mIsComputed(false)
        {
        }

        ~PluckerTreeNode()
        {
            delete mPolytope;
        }
        void setChildren(PluckerTreeNode* aLeft, PluckerTreeNode* aRight, int aSplitting)
        {
            mLeft = aLeft;
            mRight = aRight;
            mSplittingPlane = aSplitting;
        }

        PluckerPolytope<P>* getPolytope() { return mPolytope; }
        PluckerTreeNode<P>* getLeft() { return mLeft; }
        PluckerTreeNode<P>* getRight() { return mRight; }
        int getSplittingPlane()const { return mSplittingPlane; }
        bool isLeaf() const { return mLeft == nullptr && mRight == nullptr; }
        bool getIsComputed()const { return mIsComputed; }
        void setIsComputed(bool value) { mIsComputed = value; }

    private:
        PluckerPolytope<P>* mPolytope;
        int mSplittingPlane;
        PluckerTreeNode* mLeft;
        PluckerTreeNode* mRight;
        bool mIsComputed;
    };

    /** @brief Represents a complex of polytopes in Plucker space.

    The polytope complex is encoded as a BSP tree in Plucker space
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
            mRoot = new PluckerTreeNode<P>(polytope, -1);
        }

        void clear(PluckerTreeNode<P> * node);

        void computeExtremalStabbingLines(std::vector<MathVector3d> & aList, const MathPlane3d & aPlane1, const MathPlane3d & aPlane2);

        PluckerPolytope<P>* getFirstPolytope();

        void getPolytopes(std::vector<PluckerPolytope<P>*> & list)
        {
            list.clear();
            std::stack<PluckerTreeNode<P>*> s;
            s.push(mRoot);

            while (!s.empty())
            {
                PluckerTreeNode<P>* current = s.top();
                s.pop();
                if (current != nullptr)
                {
                    list.push_back(current->getPolytope());
                    s.push(current->getLeft());
                    s.push(current->getRight());
                }
            }
        }

        PluckerPolyhedron<P>* getPolyhedron() { return mPolyhedron; }
        PluckerTreeNode<P>* getRoot() { return mRoot; }
    private:

        PluckerPolyhedron<P>* mPolyhedron;
        PluckerTreeNode<P>* mRoot;
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

    template<class P>
    inline PluckerPolytope<P>* PluckerPolytopeComplex<P>::getFirstPolytope()
    {
        return mRoot->getPolytope();
    }

    template<class P>
    inline void PluckerPolytopeComplex<P>::clear(PluckerTreeNode<P>* node)
    {
        std::stack<PluckerTreeNode<P>*> s;
        s.push(node);

        while (!s.empty())
        {
            PluckerTreeNode<P>* current = s.top();
            s.pop();
            if (current != nullptr)
            {
                delete current->getPolytope();
                s.push(current->getLeft());
                s.push(current->getRight());
            }
            delete current;
        }
    }

    template<class P>
    inline void PluckerPolytopeComplex<P>::computeExtremalStabbingLines(std::vector<MathVector3d>& aList, const MathPlane3d& aPlane1, const MathPlane3d& aPlane2)
    {
        std::stack<PluckerTreeNode<P>*> s;
        s.push(mRoot);

        while (!s.empty())
        {
            PluckerTreeNode<P>* current = s.top();
            s.pop();
            if (current != nullptr)
            {
                current->getPolytope()->computeExtremalStabbingLines(aList, aPlane1, aPlane2);
                s.push(current->getLeft());
                s.push(current->getRight());
            }
        }
    }
}
