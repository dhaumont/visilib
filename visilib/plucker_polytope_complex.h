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
        
        
    PluckerPolytopeComplex(int dimension)
    {
        mDimension = dimension;
        mElements.resize(dimension - 1);
    }
    
    ~PluckerPolytopeComplex()
    {
        for (for auto myElementK: mElements)
        {
            for (auto element : myElementK)
            {
                delete element;
            }
        }
    }
    int getDimension() const
    {
        return mDimension;
    }

    void reset()
    {
        for (for auto myElementK: mElements)
        {
            for (auto element : myElementK)
            {
                element->reset();
            }
        }
    }
    
    
    const std::list<PluckerVertex*>& getVertices()
    {
        return mVertices;
    }

    void appendVertex(PluckerVertex* vertex)
    {
        mVertices.push_back(vertex);
    }

    void deleteVertex(PluckerVertex* vertex)
    {
        mVertices.remove(vertex);
   
        for (auto parent : vertex->getParents())
        {
            parent->deleteChildren(vertex);
        }
        delete vertex;
    }

    const std::list<PluckerElement*>& getElements(size_t k)
    {
        if (k == 0)
        {
            return getVertices()
        }

        return mElements[k - 1];
    }
    void appendElement(PluckerElement* element, size_t k)
    {
        if (k == 0)
        {
            appendVertex(dynamic_cast<PluckerVertex*>(element));
            return;
        }

        mElements[k - 1].push_back(element);
    }

    
    void deleteElement(PluckerElement* element, size_t k)
    {
          if (k == 0)
        {
            deleteVertex(dynamic_cast<PluckerVertex*>(element));
            return;
        }

        mElements[k - 1].remove(element);
        for (auto child : element->getChildren())
        {
            child->deleteParent(element);
        }
        for (auto parent : element->getParents())
        {
            parent->deleteChildren(element);
        }
        delete element;
    }
    private:
        int mDimension;
        std::vector<std::list<PluckerElement*>> mElements;
        set::list<PluckerVertex*> mVertices;
};



class PluckerElement
{
public:
    PluckerElement(int aRank)
    {
        mRank = aRank;   
    }
    ~PluckerElement()
    {
    }

    virtual const std::list<PluckerElement*>& getChildren() const
    {
        return mEmptyChildren;
    }
    
    virtual void appendChildren(PluckerElement* aChild)
    {
    
    }

    virtual void deleteChildren(PluckerElement* aChild)
    {
        
    }

    virtual GeometryPositionType getChildrenPosition() const
    {
        return GeometryPositionType::ON_UNKNOWN_POSITION;
    }

    
    void appendParent(PluckerElement* aParent)
    {
        mParents.push_back(aParent);
        aParent->appendChildren(this);
    }

    void deleteParent(PluckerElement* aParent)
    {
        mParents.remove(aParent);
    }

    const std::list<PluckerElement*>& getParents() const
    {
        return mParents;
    }

    
    const std::vector<size_t>& getFacets() const
    {
        return mFacets;
    }

    int getRank() const
    {
        return mRank;
    }

    private:
        
        std::list<PluckerElement*> mParents;
        std::vector<size_t> mFacets;
        int mRank;      
        bool mIntersectingQuadric;
        static const std::list<PluckerElement*> mEmptyChildren;
};

class PluckerBoundingVolume
{

};

class PluckerBoundingSphere: public PluckerBoundingVolume
{
    MathVector3d mCenter;
    double mRadius;
};

class PluckerAABB : public PluckerBoundingVolume
{
    MathVector3d mMin;
    MathVector3d mMax;
};

class PluckerInnerNode : public PluckerElement
{
    void appendChildren(PluckerElement* aChild)
    {
        mChildren.push_back(aChild);
        aChild->appendParent(this);
    }

    void deleteChildren(PluckerElement* aChild)
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
};

class PluckerVertex : public PluckerElement
{
    MathPlucker6<double>& getPlucker()
    {
        return mPlucker;
    }
    protected:
        MathPlucker6<double> mPlucker;

};

template<class P>
class PluckerEdge : public PluckerInnerNode
{
    void getVertices(PluckerVertex*& aVertex0, PluckerVertex*& aVertex1)
    {
        auto iter = mChildren.begin();
        aVertex0 = *iter;
        iter = iter++;
        aVertex1 = *iter;
    }
    
    const PluckerVertex* getVertex0() const
    {
        assert(mChildren.size() == 2);
        return static_cast<PluckerVertex*>(mChildren.front());
    }

    const PluckerVertex* getVertex1() const
    {
        assert(mChildren.size() == 2);
        return static_cast<PluckerVertex*>(mChildren.back());
    }
    const P& getVertex1() const;    
    std::vector<P> mExtremalStabbingLines;
    bool mIntersectsQudric;
};

template<class P>
class PluckerPolytopeElement : public PluckerInnerNode
{
    P mRepresentativeLine;
                               /** < @brief The ESL of the polytope, at the intersection of an edge and the Plucker Quadric*/
    std::unordered_set<Silhouette*> mSilhouettes;          /** < @brief The set of silhouettes associated to the polytope*/
};

template<class P, class S>
class PluckerInterpolatedVertex : public PluckerVertex
{
    PluckerEdge<P>* mGeneratorEdge;  
    double mAlpha, mBeta;

    void interpolate(S tolerance)
    {
        mPlucker = MathGeometry::interpolate(mAlpha,
                                             mBeta,
                                             mGeneratorEdge->getVertex0()->getPlucker(),
                                             mGeneratorEdge->getVertex1()->getPlucker(),
                                             tolerance);
    }
};

} //namespace visilib


