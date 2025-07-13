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
#include <list>
#include "math_plucker_6.h"
#include "math_geometry.h"
#include "math_combinatorial.h"
#include "geometry_position_type.h"

namespace visilib
{

    /** @brief Represents a complex of polytopes in Plucker space.

*/

    enum ELEMENT_TYPE
    {
        VERTEX = 0,
        EDGE = 1,
        FACET = 2,
        POLYTOPE = 5,
        MAXIMUM = 5
    };

    template <class P>
        class IPluckerPoint
        {
            public:
                IPluckerPoint()
                {
                    mPlucker = P::Zero();
                }
                IPluckerPoint(const IPluckerPoint& other)
                {
                    mPlucker = other.mPlucker;
                }

                P &getPlucker()
                {
                    return mPlucker;
                }
                const P &getPlucker() const
                {
                    return mPlucker;
                }

                void setPlucker(const P& a)
                {
                    mPlucker = a;
                }

            private:
                P mPlucker;
        };

    class PluckerElement
    {
        public:
            PluckerElement(int aRank)
                : mRank(aRank),
                mQuadricRelativePosition(GeometryPositionType::ON_UNKNOWN_POSITION)
            {
            }

            PluckerElement(const PluckerElement& other)
            {
                mParents = other.mParents;
                mFacetDescription = other.mFacetDescription;
                mRank = other. mRank;
                mQuadricRelativePosition = other.mQuadricRelativePosition;
            }

            ~PluckerElement()
            {
            }

            static void unlink(PluckerElement *aParent, PluckerElement *aChild)
            {
                aParent->deleteChildren(aChild);
                aChild->deleteParent(aParent);
            }

            static const std::list<PluckerElement *> &getChildren()
            {
                return mEmptyChildren;
            }

            virtual void appendChildren(PluckerElement *aChild)
            {
            }

            virtual void deleteChildren(PluckerElement *aChild)
            {
            }

            virtual GeometryPositionType getChildrenPosition() const
            {
                return GeometryPositionType::ON_UNKNOWN_POSITION;
            }

            void appendParent(PluckerElement *aParent)
            {
                mParents.push_back(aParent);
            }

            void deleteParent(PluckerElement *aParent)
            {
                mParents.remove(aParent);
            }

            const std::list<PluckerElement *> &getParents() const
            {
                return mParents;
            }

            void setFacetDescription(const std::vector<size_t> &aFacets)
            {
                mFacetDescription = aFacets;
            }

            const std::vector<size_t> &getFacetDescription() const
            {
                return mFacetDescription;
            }

            const size_t getFacetDescription(size_t index) const
            {
                return mFacetDescription[index];
            }

            GeometryPositionType getQuadricRelativePosition() const
            {
                return mQuadricRelativePosition;
            }

            void setQuadricRelativePosition(GeometryPositionType aPosition)
            {
                mQuadricRelativePosition = aPosition;
            }

            static void link(PluckerElement *aParent, PluckerElement *aChild)
            {
                aParent->appendChildren(aChild);
                aChild->appendParent(aParent);
            }


            int getRank() const
            {
                return mRank;
            }

            bool hasAncestor(PluckerElement *anAncestor) const
            {
                std::stack<const PluckerElement *> stack;
                stack.push(this);

                while (!stack.empty())
                {
                    const PluckerElement *element = stack.top();
                    stack.pop();

                    for (auto parent : element->getParents())
                    {
                        if (parent == anAncestor)
                        {
                            return true;
                        }
                        stack.push(parent);
                    }
                }
                return false;
            }

        private:
            std::list<PluckerElement *> mParents;
            std::vector<size_t> mFacetDescription;
            int mRank;
            GeometryPositionType mQuadricRelativePosition;
            static const std::list<PluckerElement *> mEmptyChildren;
    };


    class PluckerBoundingVolume
    {
        virtual bool hasIntersection() = 0;
    };

    template<class P, class S>
        class PluckerBoundingSphere: public PluckerBoundingVolume
    {
        public:
            PluckerBoundingSphere()
                : PluckerBoundingVolume()
            {

            }
            PluckerBoundingSphere(const PluckerBoundingSphere& other)
                : PluckerBoundingVolume(other)
            {
                mCenter = other.mCenter;
                mRadius = other.mRadius;
            }
            virtual ~PluckerBoundingSphere()
            {

            }
            virtual bool hasIntersection()
            {
                return false;
            }
        private:
            P mCenter;
            S mRadius;
    };

    template<class P, class S>
        class PluckerAABB : public PluckerBoundingVolume
    {
        public:
            PluckerAABB()
                : PluckerBoundingVolume()
            {

            }

            PluckerAABB(const PluckerAABB& other)
                : PluckerBoundingVolume(other)
            {
                mMin = other.mMin;
                mMax = other.mMax;
            }


            virtual ~PluckerAABB()
            {

            }
            virtual bool hasIntersection()
            {
                return false;
            }
        private:
            P mMin;
            P mMax;
    };

    class PluckerInnerNode : public PluckerElement
    {
        public:
        PluckerInnerNode(int aRank)
            : PluckerElement(aRank),
            mBoundingVolume(nullptr)
        {

        }
        PluckerInnerNode(ELEMENT_TYPE aType)
            : PluckerElement(static_cast<int>(aType)),
            mBoundingVolume(nullptr)
        {

        }

        PluckerInnerNode(const PluckerInnerNode& other)
            : PluckerElement(other)
        {
            mChildren = other.mChildren;
            mBoundingVolume = other.mBoundingVolume;
        }

        virtual void appendChildren(PluckerElement* aChild)
        {
            mChildren.push_back(aChild);
        }

        virtual void deleteChildren(PluckerElement* aChild)
        {
            mChildren.remove(aChild);
        }

        const std::list<PluckerElement*>& getChildren() const
        {
            return mChildren;
        }

        GeometryPositionType getChildrenPosition() const
        {
            GeometryPositionType position = ON_BOUNDARY;

            for (auto vertices = mChildren.begin(); vertices != mChildren.end(); vertices++)
            {

            }
            return position;
        }

        std::list<PluckerElement*> mChildren;
        PluckerBoundingVolume* mBoundingVolume;
    };

    template <class P>
        class PluckerFacet : public PluckerInnerNode, public IPluckerPoint<P>
    {
        public:
            PluckerFacet()
                : PluckerInnerNode(FACET)
            {
            }

            PluckerFacet(const PluckerFacet& other)
                : PluckerInnerNode(other),
                IPluckerPoint<P>(other)
            {
            }

            void setFacet(size_t aFacet)
            {
                std::vector<size_t> myFacets={aFacet};
                setFacetDescription(myFacets);
            }

            size_t getFacet() const
            {
                return mFacetDescription[0];
            }

    };

    template <class P>
        class PluckerVertex : public PluckerElement, public IPluckerPoint<P>
    {
        public:
        PluckerVertex()
            : PluckerElement(VERTEX)
        {
        }

        PluckerVertex(const PluckerVertex& other)
            : PluckerElement(other),
            IPluckerPoint<P>(other)
        {

        }
    };

    template <class P>
        class PluckerEdge : public PluckerInnerNode
    {
        public:
        PluckerEdge()
            : PluckerInnerNode(EDGE)
        {
        }

        PluckerEdge(const PluckerEdge& other)
            : PluckerInnerNode(other)
        {
            mExtremalStabbingLines = other.mExtremalStabbingLines;
        }

        ~PluckerEdge()
        {
        }

        void getVertices(PluckerVertex<P> *&aVertex0, PluckerVertex<P> *&aVertex1)
        {
            auto iter = mChildren.begin();
            aVertex0 = *iter;
            iter = iter++;
            aVertex1 = *iter;
        }

        const PluckerVertex<P> *getVertex0() const
        {
            assert(mChildren.size() == 2);
            return static_cast<PluckerVertex<P> *>(mChildren.front());
        }

        const PluckerVertex<P> *getVertex1() const
        {
            assert(mChildren.size() == 2);
            return static_cast<PluckerVertex<P> *>(mChildren.back());
        }

        std::vector<P> mExtremalStabbingLines;
    };

    template<class P, class S>
        class PluckerInterpolatedVertex : public PluckerVertex<P>
    {
        public:

            PluckerInterpolatedVertex(PluckerEdge<P>* aGeneratorEdge, S aAlpha, S aBeta)
                : PluckerVertex<P>(),
                mGeneratorEdge(aGeneratorEdge),
                mAlpha(aAlpha),
                mBeta(aBeta)
        {

        }
            PluckerInterpolatedVertex(const PluckerInterpolatedVertex& other)
                : PluckerVertex<P>(other)
            {
                mGeneratorEdge = other.mGeneratorEdge;
                mAlpha = other.mAlpha;
                mBeta = other.mBeta;
            }
            void interpolate(S tolerance);

        private:
            PluckerEdge<P>* mGeneratorEdge;
            S mAlpha, mBeta;

    };
} // namespace visilib
