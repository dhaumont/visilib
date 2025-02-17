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
    
    enum EntityDimension
    {
        VERTEX = 0,
        EDGE = 1,
        FACET = 4
        POLYTOPE = 5
    };

    static getFacetsCount(int dimension)
    {
        return {5, 4, 3, 2, 1, 0}[dimension];
    }    
        
    static getArrayIndex(int dimension)
    {
        return {-1, 0, 1, 2, 3, 4}[dimension];
    }    

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

    void deleteVertex(PluckerVertex* vertex)
    {
        mVertices.remove(vertex);
   
        for (auto parent : vertex->getParents())
        {
            parent->deleteChildren(vertex);
        }
        delete vertex;
    }
    
    const size_t getElementsCount(size_t k)
    {
        if (k == 0)
        {
            return getVertices().size();
        }

        return mElements[getArrayIndex(k)].size();
    }


    const std::list<PluckerElement*>& getElements(size_t k)
    {
        if (k == 0)
        {
            return getVertices()
        }

        return mElements[getArrayIndex(k)];
    }
    
    void appendVertex(PluckerVertex* vertex)
    {
        mVertices.push_back(vertex);
    }

    void appendFacet(PluckerFacet* facet)
    {        
        mElements[getArrayIndex(k)].push_back(facet);
    }

    void appendEdge(PluckerEdge* edge)
    {        
        mElements[k].push_back(edge);
    }

    void appendPolytope(PluckerPolytope* polytope)
    {
        int k = POLYTOPE - 1;
        mElements[k].push_back(polytope);
    }

    void appendElement(PluckerElement* element, size_t k)
    {
        switch(k)
        {
            case VERTEX:
                appendVertex(dynamic_cast<PluckerVertex*>(element));
                break;
            case EDGE:
                appendEdge(dynamic_cast<PluckerEdge*>(element));
                break;
            case FACET:
                appendFacet(dynamic_cast<PluckerFacet*>(element));
                break;
            case POLYTOPE:
                appendFacet(dynamic_cast<PluckerPolytope*>(element));
                break;
            default:
                mElements[k - 1].push_back(element);
        }
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

template<class P>
class IPluckerPoint
{
    public:
    P& getPlucker()
    {
        return mPlucker;
    }
    protected:
        P mPlucker;
};


class PluckerElement
{
public:
    PluckerElement(int aRank)
    : mRank(aRank)
    {
          
    }
    ~PluckerElement()
    {
    }

    static void link(PluckerElement* aParent, PluckerElement* aChild)
    {
        aParent->appendChildren(aChild);
        aChild->appendParent(aParent);
    }

    static void unlink(PluckerElement* aParent, PluckerElement* aChild)
    {
        aParent->deleteChildren(aChild);
        aChild->deleteParent(aParent);
    }

    static const std::list<PluckerElement*>& getChildren() const
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
        GeometryPositionType mPosition;
        static const std::list<PluckerElement*> mEmptyChildren;
};

class PluckerBoundingVolume
{
    virtual bool hasIntersection() = 0;
};


class PluckerBoundingSphere: public PluckerBoundingVolume
{
    virtual bool hasIntersection();
    private:
        MathVector3d mCenter;
        double mRadius;
};

class PluckerAABB : public PluckerBoundingVolume
{
    public:
        virtual bool hasIntersection();
    private:
        MathVector3d mMin;
        MathVector3d mMax;
};

class PluckerInnerNode : public PluckerElement
{
    PluckerInnerNode(int aRank)
    : PluckerElement(aRank),
    mBoundingVolume(nullptr)
    {

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



template<class P>
class PluckerFacet : public PluckerInnerNode, IPluckerPoint<P>
:: PluckerInnerNode(4)
{
    
};

template<class P>
class PluckerVertex : public PluckerElement, IPluckerPoint<P>
:: PluckerElement(0)
{
    
};

template<class P>
class PluckerEdge : public PluckerInnerNode
: PluckerInnerNode(1)
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
    bool mIntersectsQuadric;
};

template<class P>
class PluckerPolytopeElement : public PluckerInnerNode, IPluckerPoint<P>
{
    PluckerPolytopeElement()
    : PluckerInnerNode(5)
    {
        
    }
    const P& getRepresentativeLine()
    {
        return getPlucker();
    }
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

class PluckerElementFactory
{
    public:
        template<class P>
        static PluckerElement* createElement(int rank)
        {
            switch(rank)
            {
                case 0:
                    return createVertex<P>();
                case 1:
                    return createEdge<P>();
                case 4:
                    return createFacet<P>();
                case 5:
                    return createPolytope<P>();
                default:
                    return new PluckerElement<P>(rank);
            }
        }

        template<class P>
        static PluckerElement* createPolytope()
        {
            return new PluckerPolytopeElement<P>();
        }

        template<class P>
        static PluckerElement* createFacet()
        {
            return new PluckerFacet<P>();
        }

        template<class P>
        static PluckerElement* createVertex()
        {
            return new PluckerVertex<P>();
        }

        template<class P>
        static PluckerElement* createEdge()
        {
            return new PluckerEdge<P>();
        }
};

} //namespace visilib


