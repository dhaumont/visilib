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
#include "plucker_element.h"
#include "silhouette.h"

namespace visilib
{

    template<class P, class S>
        class PluckerPolytopeComplex;

    template<class P, class S>
        class PluckerPolytope : public PluckerInnerNode, public IPluckerPoint<P>
    {
        public:
            PluckerPolytope()
                : PluckerInnerNode(POLYTOPE)
            {

            }

            ~PluckerPolytope()
            {
                // TODO: delete mSilhouette<P>s
            }

            const P& getRepresentativeLine() const
            {
                return this->getPlucker();
            }

            const PluckerPolytopeComplex<P,S>* getPolytopeComplex() const
            {
                assert(mParents.size() == 1);
                return static_cast<PluckerPolytopeComplex<P,S>*>(mParents.front());
            }

            void fillEdges(std::vector<PluckerEdge<P>*>& aPolytopeEdges);
            void fillVertices(std::vector<PluckerVertex<P>*>& aPolytopeVertices);
            bool containsRealLines() const;
            void getExtremalStabbingLinesBackTo3D(std::vector<std::pair<MathVector3d, MathVector3d>>& aStabbingLines, const MathPlane3d& aPlane1, const MathPlane3d& aPlane2);
            void computeEdgesIntersectingQuadric(S tolerance);
        private:
            /** < @brief The ESL of the polytope, at the intersection of an edge and the Plucker Quadric*/
            std::unordered_set<Silhouette<P>*> mSilhouettes;          /** < @brief The set of silhouettes associated to the polytope*/
    };



    template<class P, class S>
        class PluckerPolytopeComplex : public PluckerInnerNode
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

            int getCombinatorialFacetsMaximumCount(int dimension) const
            {
                return getDimension() - dimension;
            } 

            int getDimension() const
            {
                return mElements.size();
            }

            const size_t getElementsCount(size_t k)
            {
                return mElements[k].size();
            }

            const std::list<PluckerElement*>& getElements(size_t k)
            {
                return mElements[k];
            }

            void appendVertex(PluckerVertex<P>* vertex)
            {
                mElements[VERTEX].push_back(vertex);
            }

            void appendFacet(PluckerFacet<P>* facet)
            { 
                mElements[FACET].push_back(facet);
            }

            void appendEdge(PluckerEdge<P>* edge)
            { 
                mElements[EDGE].push_back(edge);
            }

            void appendPolytope(PluckerPolytope<P,S>* polytope)
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

            void deleteVertex(PluckerVertex<P>* vertex)
            {
                deleteElement(vertex,VERTEX);
            }

            void deleteEdge(PluckerEdge<P>* edge)
            {
                deleteElement(edge,EDGE);
            }

            void deleteFacet(PluckerFacet<P>* facet)
            {
                deleteElement(facet,FACET);
            }

            void deletePolytope(PluckerPolytope<P,S>* polytope)
            {
                deleteElement(polytope,POLYTOPE);
            }

            template<class T>
                class ElementIterator
                {
                    public:
                        ElementIterator(const std::list<PluckerElement*>::iterator& aIterator, PluckerElement* anAncestor = nullptr)
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
                    private:
                        PluckerElement* mAncestor;
                        typename std::list<T*>::iterator mIterator;
                };

            typedef ElementIterator<PluckerVertex<P>> VertexIterator;
            typedef ElementIterator<PluckerEdge<P>> EdgeIterator;
            typedef ElementIterator<PluckerFacet<P>> FacetIterator;
            typedef ElementIterator<PluckerPolytope<P,S>> PolytopeIterator;


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
            template<class P, class S>
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
                            return new PluckerElement(rank);
                    }
                }

            template<class P, class S>
                static PluckerElement* createPolytope()
                {
                    return new PluckerPolytope<P,S>();
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

    template<class P,class S>
        void PluckerPolytope<P,S>::fillVertices(std::vector<PluckerVertex<P>*>& aPolytopeVertices)
        {
            PluckerPolytopeComplex<P,S>* complex = getPolytopeComplex();

            for (auto iter = complex->beginVertices(this); iter != complex->endVertices(this); iter++)
            {
                aPolytopeVertices.push_back(*iter);
            } 
        }

    template<class P, class S>
        void PluckerPolytope<P,S>::fillEdges(std::vector<PluckerEdge<P>*>& aPolytopeEdges)
        {
            PluckerPolytopeComplex<P,S>* complex = getPolytopeComplex();

            for (auto iter = complex->beginEdge(this); iter != complex->endEdge(this); iter++)
            {
                aPolytopeEdges.push_back(*iter);
            } 
        }

    template<class P, class S>
        bool PluckerPolytope<P,S>::containsRealLines() const
        {
            PluckerPolytopeComplex<P,S>* complex = getPolytopeComplex();

            int myCount = 0;
            for (auto iter = complex->beginEdges(this); iter != complex->endEdges(this); iter++)
            {
                PluckerEdge<P>* myEdge = *iter;

                if (myEdge->getQuadricRelativePosition() == ON_BOUNDARY)
                {
                    return true;
                }
                myCount = myCount++;
            }

            return myCount > 0;
        } 

    template<class P, class S>
        void PluckerPolytope<P,S>::getExtremalStabbingLinesBackTo3D(std::vector<std::pair<MathVector3d, MathVector3d>>& aStabbingLines, const MathPlane3d& aPlane1, const MathPlane3d& aPlane2)
        {
            aStabbingLines.clear();

            PluckerPolytopeComplex<P,S>* complex = getPolytopeComplex();
            std::vector<PluckerEdge<P>*> myEdges;
            fillEdges(myEdges);

            for (auto iter : myEdges)
            {
                PluckerEdge<P>* myEdge = *iter;

                if (myEdge->getQuadricRelativePosition() == ON_BOUNDARY)
                {
                    std::pair<MathVector3d, MathVector3d> r = MathGeometry::getBackTo3D(myEdge->getPluckerPoint(), aPlane1, aPlane2);
                    aStabbingLines.push_back(r);
                }
            }

#if 1

            for (size_t i = 0; i < myEdges.size(); i++)
            {
                for (size_t j = i + 1; j < myEdges.size(); j++)
                {
                    if (MathCombinatorial::haveAtLeastNCommonFacets(myEdges[i]->getFacetDescription(), myEdges[j]->getFacetDescription(), 3))
                    {
                        aStabbingLines.push_back(std::pair<MathVector3d, MathVector3d>( aStabbingLines[i].first,  aStabbingLines[j].first));
                        aStabbingLines.push_back(std::pair<MathVector3d, MathVector3d>( aStabbingLines[i].second, aStabbingLines[j].second));
                    }
                }
            }
#endif

            V_ASSERT(aStabbingLines.size() > 0);
        }

    template<class P,class S>
        void PluckerPolytope<P,S>::computeEdgesIntersectingQuadric(S tolerance)
        {
            std::vector<PluckerEdge<P>*> myEdges;
            fillEdges(myEdges);

            if (myEdges.size() > 0 && MathPredicates::hasPluckerPolytopeIntersectionWithQuadric(this))
            {
                for (auto iter : myEdges)
                {
                    PluckerEdge<P>* myEdge = *iter; 
                    if (myEdge->getQuadricRelativePosition() == UNKNOWN)
                    {
                        //At least one
                        PluckerVertex<P>* vertex0 = myEdge->getVertex0(); 
                        PluckerVertex<P>* vertex1 = myEdge->getVertex1(); 
                        {
                            const P& v0 = vertex0->getPlucker();
                            const P& v1 = vertex1->getPlucker();
                            GeometryPositionType p0 = vertex0->getQuadricRelativePosition();
                            GeometryPositionType p1 = vertex1->getQuadricRelativePosition();
                            if (MathGeometry::hasPluckerEdgeWithQuadricIntersection(v0, v1, p0, p1, tolerance))
                            {
                                myEdge->setQuadricRelativePosition(ON_BOUNDARY);
                            }
                        }
                    }
                }
            }
        }
} //namespace visilib

