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
#include "geometry_position_type.h"
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

   enum ELEMENT_TYPE
   {
       VERTEX = 0,
       EDGE = 1,
       FACET = 2,
       POLYTOPE = 5,
       MAXIMUM = 5
   };

       
template<class P>
class IPluckerPoint
{
    public:
    P& getPlucker()
    {
        return mPlucker;
    }
    const P& getPlucker() const
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
    : mRank(aRank),
      mQuadricRelativePosition(GeometryPositionType::ON_UNKNOWN_POSITION)
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

    static const std::list<PluckerElement*>& getChildren()
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

    void setFacetDescription(const std::vector<size_t>& aFacets)
    {
        mFacetDescription = aFacets;
    }

    const std::vector<size_t>& getFacetDescription() const
    {
        return mFacetDescription;
    }

    const size_t setFacetDescription(size_t index) const
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

    int getRank() const
    {
        return mRank;
    }

    bool hasAncestor(PluckerElement* anAncestor) const
    {
        std::stack<const PluckerElement*> stack;
        stack.push(this);
         
        while (!stack.empty())
        {
            const PluckerElement* element = stack.top();
            stack.pop();
            
            for (auto parent: element->getParents())
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
        
        std::list<PluckerElement*> mParents;
        std::vector<size_t> mFacetDescription;
        int mRank;      
        GeometryPositionType mQuadricRelativePosition;
        static const std::list<PluckerElement*> mEmptyChildren;
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
:: PluckerInnerNode(FACET)
{
    
};

template<class P>
class PluckerVertex : public PluckerElement, IPluckerPoint<P>
:: PluckerElement(VERTEX)
{
    
};

template<class P>
class PluckerEdge : public PluckerInnerNode
{
    PluckerEdge()
    : PluckerInnerNode(EDGE)
    {

    }

    ~PluckerEdge()
    {

    }

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
};

template<class P>
class PluckerPolytopeComplex;

template<class P>
class PluckerPolytope : public PluckerInnerNode, IPluckerPoint<P>
{
    PluckerPolytope()
    : PluckerInnerNode(POLYTOPE)
    {
        
    }

    ~PluckerPolytope()
    {
        // TODO: delete mSilhouettes
    }
    const P& getRepresentativeLine()
    {
        return getPlucker();
    }
    
    void fillVertices(PluckerPolytopeComplex<P>* complex, std::vector<PluckerVertex<P>*>& aPolytopeVertices);
    


                               /** < @brief The ESL of the polytope, at the intersection of an edge and the Plucker Quadric*/
    std::unordered_set<Silhouette*> mSilhouettes;          /** < @brief The set of silhouettes associated to the polytope*/
};

template<class P, class S>
class PluckerInterpolatedVertex : public PluckerVertex
{
    public:

    PluckerInterpolatedVertex(PluckerEdge<P>* aGeneratorEdge, double aAlpha, double aBeta)
    : PluckerVertex(),
      mGeneratorEdge(aGeneratorEdge),
      mAlpha(aAlpha),
      mBeta(aBeta)
    {
     
    }
    template<class S>
    void interpolate(S tolerance)
    {
        mPlucker = MathGeometry::interpolate(mAlpha,
                                             mBeta,
                                             mGeneratorEdge->getVertex0()->getPlucker(),
                                             mGeneratorEdge->getVertex1()->getPlucker(),
                                             tolerance);
    }
    private:
        PluckerEdge<P>* mGeneratorEdge;  
        double mAlpha, mBeta;

};
 
    template<class P>
    class PluckerPolytopeComplex
    {
    public:
    

    PluckerPolytopeComplex(int dimension = MAXIMUM)    
    {
        mElements.resize(dimension);
    }

    ~PluckerPolytopeComplex()
    {
        for (auto myElementK: mElements)
        {
            for (auto element : myElementK)
            {
                delete element;
            }
        }
    }

    static getCombinatorialFacetsMaximumCount(int dimension)
    {
        return {5, 4, 3, 2, 1, 0}[dimension];
    }    
    
    int getDimension() const
    {
        return mElements.size();
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
        
    const size_t getElementsCount(size_t k)
    {
        return mElements[k].size();
    }

    const std::list<PluckerElement*>& getElements(size_t k)
    {
        return mElements[k];
    }
    
    void appendVertex(PluckerVertex* vertex)
    {
        mElements[VERTEX].push_back(vertex);
    }

    void appendFacet(PluckerFacet* facet)
    {        
        mElements[FACET].push_back(facet);
    }

    void appendEdge(PluckerEdge* edge)
    {        
        mElements[EDGE].push_back(edge);
    }

    void appendPolytope(PluckerPolytope* polytope)
    {        
        mElements[POLYTOPE].push_back(polytope);
    }

    void appendElement(PluckerElement* element, size_t k)
    {
        mElements[k].push_back(element);       
    }

    void deleteElement(PluckerElement* element, size_t k)
    {           
        for (auto child : element->getChildren())
        {
            child->deleteParent(element);
        }
        for (auto parent : element->getParents())
        {
            parent->deleteChildren(element);
        }
        mElements[k].remove(element);
        delete element;
    }

    void deleteVertex(PluckerVertex* vertex)
    {
        deleteElement(vertex,VERTEX);
    }

    void deleteEdge(PluckerEdge* edge)
    {
        deleteElement(edge,EDGE);
    }

    void deleteFacet(PluckerFacet* facet)
    {
        deleteElement(facet,FACET);
    }

    void deletePolytope(PluckerPolytope* polytope)
    {
        deleteElement(polytope,POLYTOPE);
    }

    template<class T>
    struct ElementIterator
    {
        ElementIterator(std::list<PluckerElement*>::iterator aIterator, PluckerElement* anAncestor = nullptr)
        :  mIterator(aIterator),
           mAncestor(anAncestor)            
        {            
        }

        T operator*()
        {
            return static_cast<T*>(*mIterator);
        }

        void operator++()
        {
            if (mAncestor == nullptr)
            {
                mIterator++;
            }
            else
            {
                while (mIterator != mElements.end() && (*mIterator)->getAncestor() != mAncestor)
                {
                    mIterator++;
                }
            }
        }

        bool operator!=(const T& aIterator)
        {
            return mIterator != aIterator.mIterator;
        }
        
        std::list<T*>::iterator mIterator;
    };

    typedef ElementIterator<PluckerVertex> VertexIterator;
    typedef ElementIterator<PluckerEdge> EdgeIterator;
    typedef ElementIterator<PluckerFacet> FacetIterator;
    typedef ElementIterator<PluckerPolytope> PolytopeIterator;

    
    VertexIterator beginVertices(PluckerElement* anAncestor = nullptr)
    {
        return VertexIterator(mElements[VERTEX].begin(), anAncestor);
    }
    VertexIterator endVertices(PluckerElement* anAncestor = nullptr)
    {
        return VertexIterator(mElements[VERTEX].end(), anAncestor);
    }

    EdgeIterator beginEdges(PluckerElement* anAncestor = nullptr)
    {
        return EdgeIterator(mElements[EDGE].begin(), anAncestor);
    }

    EdgeIterator endEdges(PluckerElement* anAncestor = nullptr)
    {
        return EdgeIterator(mElements[EDGE].end(), anAncestor);
    }

    FacetIterator beginFacets(PluckerElement* anAncestor = nullptr)
    {
        return FacetIterator(mElements[FACET].begin(), anAncestor);
    }

    FacetIterator endFacets(PluckerElement* anAncestor = nullptr)
    {
        return FacetIterator(mElements[FACET].end(), anAncestor);
    }

    PolytopeIterator beginPolytopes(PluckerElement* anAncestor = nullptr)
    {
        return PolytopeIterator(mElements[POLYTOPE].begin(), anAncestor);
    }

    PolytopeIterator endPolytopes(PluckerElement* anAncestor = nullptr)
    {
        return PolytopeIterator(mElements[POLYTOPE].end(), anAncestor);
    }
    private:        
        std::vector<std::list<PluckerElement*>> mElements;
};

class PluckerElementFactory
{
    public:
        template<class P>
        static PluckerElement* createElement(int rank)
        {
            switch(rank)
            {
                case VERTEX:
                    return createVertex<P>();
                case EDGE:
                    return createEdge<P>();
                case FACET:
                    return createFacet<P>();
                case POLYTOPE:
                    return createPolytope<P>();
                default:
                    return new PluckerElement<P>(rank);
            }
        }

        template<class P>
        static PluckerElement* createPolytope()
        {
            return new PluckerPolytope<P>();
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

template<class P>
void PluckerPolytope<P>::fillVertices(PluckerPolytopeComplex<P>* complex, std::vector<PluckerVertex<P>*>& aPolytopeVertices)
{
    for (auto iter = complex->beginVertices(this); iter != complex->endVertices(this); iter++)
    {
        aPolytopeVertices.push_back(*iter);
    }    

} //namespace visilib


