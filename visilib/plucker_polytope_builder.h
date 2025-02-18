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

#include "math_predicates.h"
#include "geometry_convex_polygon.h"
#include "math_combinatorial.h"
#include "plucker_polytope.h"

namespace visilib
{
    /** @brief Create the minimal polytope containing all the stabbing lines for polygon a and b */

    template <class P, class S>
    class PluckerPolytopeBuilder
    {
    public:
        PluckerPolytopeBuilder(bool normalization, S tolerance);

        /** @brief Create the minimal polytope following Mora et al. An extermal stabbing line is created for each vertex pair of polygon a and b.
        @param a: first source polygon
        @param b: second source polygon
        @param polyhedron: the polyhedron that will contains all the vertices of the polytope
        @return: the created polytope
        */
        PluckerPolytopeComplex<P>* build(GeometryConvexPolygon& a, GeometryConvexPolygon& b);

    private:
        /** @brief Add the edges of polygon a and b as hyperplanes of the polyhedron. The approximate normal is used when the polygon is degenerated */
        void addSourcePolygonEdgesAsHyperplanes(PluckerPolytopeComplex<P>* polytope, const GeometryConvexPolygon& aPolygon, bool clockWise, const MathVector3d& approximateNormal);

        /** @brief Add the extremal stabbing lines of a and b as hyperplanes of the polyhedron */
        void addExtremalStabbingLines(PluckerPolytopeComplex<P>* polytope, GeometryConvexPolygon& a, GeometryConvexPolygon& b);

        /** @brief Create the edges from the facets representation to obtain the polytope skeleton*/
        void addEdges(PluckerPolytopeComplex<P>* polytope);

        void addStabbingLinesToEdges(PluckerPolytopeComplex<P> * polytope);

        bool mNormalization;  /** < @brief Normalization activated or not*/
        S mTolerance;         /** < @brief Tolerance used to compute the intersection with the Plucker quadric */
    };

    template<class P, class S>
    inline PluckerPolytopeBuilder<P, S>::PluckerPolytopeBuilder(bool aNormalization, S aTolerance)
    {
        mNormalization = aNormalization;
        mTolerance = aTolerance;
    }

    template<class P, class S>
    inline void PluckerPolytopeBuilder<P, S>::addSourcePolygonEdgesAsHyperplanes(PluckerPolytopeComplex<P>* polytope, const GeometryConvexPolygon& aPolygon, bool clockWise, const MathVector3d& approximateNormal)
    {
        size_t n = aPolygon.getVertexCount();
        std::vector<MathVector3d> v;

        if (n < 3)
        {
            // manage the case of a degenerated polygons
            v.resize(3);
            v[1] = aPolygon.getVertex(0);

            if (n == 1)
            {
                // the polygon is a single vertex. We initialize two temporary vertices
                v[0] = v[1] + MathVector3d::X();
                v[2] = v[1] + MathVector3d::Y();
            }
            else if (n == 2)
            {
                // the polygon is a single edge. We initialize one temporary vertex
                v[2] = aPolygon.getVertex(1);

                MathVector3d edge = v[2] - v[1];
                MathVector3d e1 = MathVector3d::cross(approximateNormal, edge);
                e1.normalize();
                e1 *= edge.getNorm();

                v[0] = (v[1] + v[2]) * 0.5 + e1;
            }

            /* We compute an additional vertex that will be used to add two edges to form a valid polygon, in order to make sure that each vertex of the polygon
                is at the intersection of two edges. This ensure later one in Plucker space that each extremal lines is at the intersection of four hyperplanes.*/

            MathVector3d edge = v[2] - v[1];
            MathVector3d e1 = MathVector3d::cross(approximateNormal, edge);
            e1.normalize();
            e1 *= edge.getNorm();

            v[0] = (v[1] + v[2]) * 0.5 + e1;
        }
        else
        {
            // manage the case of a regular convex polygon
            v = aPolygon.getVertices();
        }
        size_t count = v.size();

        // Compute the oriented Plucker line hyperplane representation of each polygonal edge and store it
        for (size_t i = 0; i < count; i++)
        {
            const MathVector3d& v1 = v[i % count];
            const MathVector3d& v2 = v[(i + 1) % count];

            P  line = clockWise ? P(v1, v2) : P(v2, v1);
            if (mNormalization)
            {
                line = line.getNormalized();

                V_ASSERT((MathPredicates::isNormalized(line, mTolerance)));
            }
            PluckerFacet<P>* myFacet = new PluckerFacet<P>(line);
            myFacet.mFacets.push_back(polytope->getElementsCount(4)}
            polytope->appendFacet(myFacet);
            
    }

    template<class P, class S>
    inline void PluckerPolytopeBuilder<P, S>::addExtremalStabbingLines(PluckerPolytopeComplex<P> * polytope, GeometryConvexPolygon & a, GeometryConvexPolygon & b)
    {
        // There is an extremal stabbing line passing through each vertex of polygon a and each vertex of polygon b.
        // This procedure initialize the facets of each stabbing line, which is the list of edges incident to the lines in 3D space. In Plucker space, it is list of hyperplanes that meet at
        // the Plucker point corresponding to the stabbing line
                
        size_t count1 = a.getVertexCount();
        size_t count2 = b.getVertexCount();

        size_t mySize1 = count1 > 3 ? count1 : 3;

        std::vector<size_t> myTemp(4);

        int myStabbingNumber = 0;

        // we iterate through all the vertex of polygon a
        for (size_t i = 0; i < count1; i++)
        {
            if (count1 > 2)
            {
                // regular case: the first two facets are the edges of the polygon a
                MathVector2i edges = a.getEdgesOfVertex(i);
                myTemp[0] = edges.x;
                myTemp[1] = edges.y;
            }
            else
            {
                // degenated case (polygon a is a single vertex or an edge): the facets are the additional edges that we have added
                myTemp[0] = i;
                myTemp[1] = i + 1;
            }

            // we iterate through all the vertex of polygon b
            for (int j = 0; j < count2; j++)
            {
                const MathVector3d& v1 = a.getVertex(i);
                const MathVector3d& v2 = b.getVertex(j);
                P  myStabbingLine(v1, v2);
                if (mNormalization)
                {
                    myStabbingLine = myStabbingLine.getNormalized();
                }
                
                PluckerVertex<P>* myVertex = new PluckerVertex<P>(myStabbingLine);
                myVertex->setQuadricRelativePosition(ON_BOUNDARY);
                polytope->appendVertex(myVertex);

                if (count2 > 2)
                {
                    // regular case: the first two facets are the edges of the polygon b
                    MathVector2i edges = b.getEdgesOfVertex(j);

                    myTemp[2] = edges.x + mySize1;
                    myTemp[3] = edges.y + mySize1;
                }
                else
                {
                    // degenated case (polygon b is a single vertex or an edge): the facets are the additional edges that we have added
                    myTemp[2] = j + mySize1;
                    myTemp[3] = j + mySize1 + 1;
                }
                myEdge->setFacetDescription(myTemp);
                
                myStabbingNumber++;
            }
        }
    }

    template<class P, class S>
    inline PluckerPolytopeComplex<P>* PluckerPolytopeBuilder<P, S>::build(GeometryConvexPolygon & a, GeometryConvexPolygon & b)
    {
        PluckerPolytopeComplex<P>* myPolytope = new PluckerPolytopeComplex<P>();

        // The gravity centers are used to compute an approximated normal vector, that will be used when polygons are degenerated
        MathVector3d ga = MathGeometry::getGravityCenter(a);
        MathVector3d gb = MathGeometry::getGravityCenter(b);

        addSourcePolygonEdgesAsHyperplanes(myPolytope, a, true, gb - ga);
        addSourcePolygonEdgesAsHyperplanes(myPolytope, b, false, ga - gb);

        addExtremalStabbingLines(myPolytope, a, b);
        addEdges(myPolytope);


        myPolytope->computeEdgesIntersectingQuadric(mTolerance);
        V_ASSERT(myPolytope->isValid(aPolyhedron, mNormalization, mTolerance));

        return myPolytope;
    }

    template<class P, class S>
    inline void PluckerPolytopeBuilder<P, S>::addEdges(PluckerPolytopeComplex<P> * polytope)
    {
        // We iterate through all the extremal stabbing lines, and creates and edge for the vertices that share at least three common facets
        // The complexity is O(n^2), and more optimized version could be devised but we keep this one for simplicity
        
        std::vector<PluckerFacet*> myFacetsVector;
        
        PluckerPolytope* myPolytope = new PluckerPolytope();
        polytope->appendPolytope(myPolytope);

        for (auto iter = polytope->beginFacets(); iter != polytope->endFacets(); iter++)
        {
            myFacetsVector.push_back(*iter);
            PluckerElement::link(myPolytope, *iter);
        }

        for (int k = 1; k < polytope->getDimension(); k++)
        {            
            for (auto iter1 : polytope->getElements(k-1))
            {
                PluckerElement* child1 = *iter1;

                for (auto iter2 = iter1; iter2 != children.end(); iter2++)
                {   
                    PluckerElement* child2 = *iter2;
                    V_ASSERT (child1!=child2)
                    
                    std::vector<size_t> myTemp;
                    MathCombinatorial::getCommonFacets(child1->getFacetDescription(), child2->getFacetDescription(), myTemp);
                    if (myTemp.size() >= getCombinatorialFacetsMaximumCount(k))
                    {
                        PluckerElement* parent;                            

                        if (k == FACET)
                        {
                            V_ASSERT(myTemp.size() == 1);

                            parent = myFacetsVector[myTemp[0]];
                            
                            V_ASSERT(myTemp[0] == parent->getFacetDescription(0));
                        }
                        else
                        {
                            parent = PluckerElementFactory::createElement(k);
                            parent->setFacetDescription(myTemp);
                            parent->setQuadricRelativePosition(ON_BOUNDARY);
                        }
                        
                        PluckerElement::link(parent, child1);
                        PluckerElement::link(parent, child2);
                        myLattice->appendElement(parent,k);
                    }
                
                }
            }
        }
    }

    template<class P, class S>
    inline void PluckerPolytopeBuilder<P, S>::addStabbingLinesToEdges(PluckerPolytopeComplex<P> * polytope)
    {        
        for (auto iter = polytope->beginEdges(); iter != polytope->endEdges(); iter++)
        {
            PluckerEdge<P>* myEdge = *iter;
            myEdge->mExtremalStabbingLines.push_back(myEdge->getVertex0()->getPlucker());
            myEdge->mExtremalStabbingLines.push_back(myEdge->getVertex1()->getPlucker());
        }
    }
}
