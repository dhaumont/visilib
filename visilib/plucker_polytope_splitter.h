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

#include <unordered_map>
#include <list>
#include <stack>
#include <map>
#include "geometry_position_type.h"
#include "math_predicates.h"
#include "math_geometry.h"
#include "math_combinatorial.h"
#include "plucker_polytope.h"

namespace visilib
{
    /** @brief Split a convex polytope into two polytopes, following an hyperplane. 
    
    The polytopes are represented by their skeleton (myVertices and edges) only. The splitting is based on the combinatorial facet description of the vertices and edges */

    template<class P, class S>
    class PluckerPolytopeSplitter
    {
    public:
        /** @brief Split a polytope by an hyperplane
        @param polyhedron: the polyhedron containing the vertices of the polytopes
        @param aPlane: the splitting hyperplane equation
        @param aPolytope: the input polytope
        @param aPlaneId: the hyperplane facet number (index of the hyperplane in the polyhedron)
        @param aNormalization: boolean to perform a normalization or not
        @param tolerance: a tolerance to determine  wheter a point lies on the splitting hyperplane or not
        @param aLeft: the resulting splitted polytope at the negative side of the hyperplane (must be instanciated before calling the function)
        @param aRight: the resulting splitted polytope at the positive side of the hyperplane  (must be instanciated before calling the function)
       */
        static GeometryPositionType split(PluckerPolyhedron<P>* polyhedron, const P& aPlane, PluckerPolytope<P>* aPolytope, PluckerPolytope<P>* aLeft, PluckerPolytope<P>* aRight, size_t aPlaneId, bool anormalization, S tolerance);
    };

    template<class P, class S>
    inline GeometryPositionType PluckerPolytopeSplitter<P, S>::split(PluckerPolyhedron<P>* aPolyhedron, const P& aPlane, PluckerPolytope<P>* aPolytope, PluckerPolytope<P>* aLeft, PluckerPolytope<P>* aRight, size_t aPlaneID, bool normalization, S tolerance)
    {
        bool hasPointOnTheLeft = false;
        bool hasPointOnTheRight = false;

        std::vector<size_t> myWait;
        std::unordered_set<size_t> myVertices;
        std::vector<size_t> myQueryList;

        std::unordered_map<size_t, int> myStatus;
        std::unordered_map<size_t, S>	myArrayOfPlanePosition;

        V_ASSERT(aPolytope->getEdgeCount() > 0);

        size_t edge[2];

         
        // Step 1 - edge Position

        // We iterate through all the edges to inspect the relative position of the two myVertices of the edge with regards to the splitting hyperplane

        for (auto iter = aPolytope->getEdges().begin(); iter != aPolytope->getEdges().end(); iter++)
        {
            edge[0] = iter->first;
            edge[1] = iter->second;

            for (size_t i = 0; i < 2; i++)
            {
                size_t a = edge[i];

                if (myVertices.find(a) == myVertices.end())
                {
                    myVertices.insert(a);

                    S result = aPlane.dot(aPolyhedron->get(a));
                    GeometryPositionType position = MathPredicates::getRelativePosition(S, tolerance);

                    myArrayOfPlanePosition.insert(std::pair<size_t, S>(a, result));

                    if (position == ON_NEGATIVE_SIDE)
                    {
                        myStatus.insert(std::pair<size_t, int>(a, -1));
                        hasPointOnTheLeft = true;
                    }
                    else if (position == ON_POSITIVE_SIDE)
                    {
                        myStatus.insert(std::pair<size_t, int>(a, 1));
                        hasPointOnTheRight = true;
                    }
                    else
                    {
                        V_ASSERT(position == ON_BOUNDARY);
                        myStatus.insert(std::pair<size_t, int>(a, 0));
                        myWait.push_back(a);
                        myQueryList.push_back(a);
                    }
                }
            }
        }

        // Step 2 - Treat the special cases, when the hyperplane does not split with the polytope

        if (!hasPointOnTheLeft && !hasPointOnTheRight)
        {
            // Step 2 - Case 1 - the hyperplanes contains all the myVertices of the polytope
            //			   	   - no split is required, the hyperplane is added to the facets description of the myVertices contained in the hyperplane.

            for (auto iter = aPolytope->getEdges().begin(); iter != aPolytope->getEdges().end(); iter++)
            {
                size_t myI1 = iter->first;
                size_t myI2 = iter->second;
                aLeft->addEdge(myI1, myI2, aPolyhedron);
                aRight->addEdge(myI1, myI2, aPolyhedron);
            }
            V_ASSERT(aLeft->getEdgeCount() > 0);
            V_ASSERT(aRight->getEdgeCount() > 0);
            for (size_t index = 0; index < myWait.size(); index++)
            {
                if (!MathCombinatorial::hasFacet(aPolyhedron->getFacetsDescription(myWait[index]), aPlaneID))
                {
                    aPolyhedron->addToFacetsDescription(myWait[index], aPlaneID);
                }
            }
            return ON_BOUNDARY;
        }
        if (!hasPointOnTheLeft)
        {
            // Step 2 - Case 2  - the polytope is entierly on the positive side. Nothing to do.
            return ON_POSITIVE_SIDE;
        }
        else if (!hasPointOnTheRight)
        {
            // Step 2 - Case 3  - the polytope is on entierly the negative side. Nothing to do.
            return ON_NEGATIVE_SIDE;
        }

        
        // Step 3 - P processing
        //			We are in the generic case, where the hyperplane split the polytope in two convex polytopes.
        //          We first process all the myVertices and creates the new myVertices if required
        //

        // The hyperplane is added to the facets description of the myVertices contained in the hyperplane.
        for (size_t index = 0; index < myWait.size(); index++)
        {
            if (!MathCombinatorial::hasFacet(aPolyhedron->getFacetsDescription(myWait[index]), aPlaneID))
            {
                aPolyhedron->addToFacetsDescription(myWait[index], aPlaneID);
            }
        }

        // We iterate though all the edges to determine the operation required for each edge
        size_t initialNumberOfVertices = aPolyhedron->getLinesCount();

        for (auto iter = aPolytope->getEdges().begin(); iter != aPolytope->getEdges().end(); iter++)
        {
            size_t myI1 = iter->first;
            size_t myI2 = iter->second;

            int sum = myStatus[myI1] + myStatus[myI2];

            //  myStatus[myI1]   |   + + 0      | + 0 - |  0  -  -    |  position of vertex 0
            //  myStatus[myI2]   |   + 0 +      | - 0 + |  -  0  -    |  position of vertex 1
            //                  |   =======    |=======|==========   |
            // sum              |   2 1 1      | 0 0 0 | -1 -1 -2    |  the sum
            // operation	    | add to right | split | add to left |  the operation to perform

            if (sum < 0) // add edge to left polytope
            {
                aLeft->addEdge(myI1, myI2, aPolyhedron);
            }
            else if (sum > 0)// add edge to right polytope
            {
                aRight->addEdge(myI1, myI2, aPolyhedron);
            }
            else // split edge
            {
                if (myStatus[myI1] != 0)
                {
                    // the edge crosses the hyperplane: a new vertex is created at the intersection the edge and the hyperplane

                    if (myStatus[myI1] > 0)
                    {
                        size_t temp = myI1;
                        myI1 = myI2;
                        myI2 = temp;
                    }
                    
                   const std::vector<size_t>& facetsI1 = aPolyhedron->getFacetsDescription(myI1);
                   const std::vector<size_t>& facetsI2 = aPolyhedron->getFacetsDescription(myI2);
         
                    std::vector<size_t> myFacets;
                    MathCombinatorial::initFacets(facetsI1, 
                                                  facetsI2, 
                                                  aPlaneID,
                                                  myFacets);


                    bool hasVertex = false;
                    for (size_t test = initialNumberOfVertices; test < aPolyhedron->getLinesCount() && !hasVertex; test++)
                    {
                        const std::vector<size_t> existingFacets = aPolyhedron->getFacetsDescription(test);
                        if (existingFacets.size() == myFacets.size() && std::equal(existingFacets.begin(), existingFacets.end(), myFacets.begin()))
                            hasVertex = true;
                    }

                    if (hasVertex)
                        continue;

                    // Compute vertex position, at the intersection of the edge and the splitting hyperplane
                    P  myIntersection = MathGeometry::interpolate(myArrayOfPlanePosition[myI1], myArrayOfPlanePosition[myI2], aPolyhedron->get(myI1), aPolyhedron->get(myI2), tolerance);
                    if (normalization)
                    {
                        myIntersection = myIntersection.getNormalized();
                    }

                    size_t vertexIndex = aPolyhedron->add(myIntersection, MathPredicates::getQuadricRelativePosition(myIntersection, tolerance), normalization, tolerance);

                    // Prepare for new edges creation
                    myQueryList.push_back(vertexIndex);

                    aPolyhedron->initFacetsDescription(vertexIndex, myFacets);

                    V_ASSERT(aPolyhedron->checkFacetsDescription(vertexIndex, myI1, myI2, aPlaneID));

                    // Add the edges to the left and to the right
                    aLeft->addEdge(myI1, vertexIndex, aPolyhedron);
                    aRight->addEdge(vertexIndex, myI2, aPolyhedron);
                }
                else
                {
                    // the myVertices are on the hyperplane, we do not have to split the edge

                    aLeft->addEdge(myI1, myI2, aPolyhedron);
                    aRight->addEdge(myI1, myI2, aPolyhedron);
                }
            }
        }

        //** @brief
        // Step 4 - Creation of the edges of the polytope skeleton

        // We iterate through all the new extremal stabbing lines, and creates and edge for the myVertices that share at least three common facets and if the edge will not be
        // degenerated (length 0)
        for (size_t m = 0; m < myQueryList.size(); m++)
        {
            for (size_t n = m + 1; n < myQueryList.size(); n++)
            {
                size_t Qm = myQueryList[m];
                size_t Qn = myQueryList[n];

                const std::vector<size_t>& facetsQm = aPolyhedron->getFacetsDescription(Qm);
                const std::vector<size_t>& facetsQn = aPolyhedron->getFacetsDescription(Qn);
         
                if (MathCombinatorial::haveAtLeastNCommonFacets(facetsQm, facetsQn))
                {
                    V_ASSERT(Qm != Qn);
                    if (!MathPredicates::isEdgeCollapsed(aPolyhedron->get(Qn), aPolyhedron->get(Qm), tolerance))
                    {
                        if (facetsQm.size() != facetsQn.size() || !std::equal(facetsQm.begin(), facetsQm.end(), facetsQn.begin()))
                        {
                        aLeft->addEdge(Qm, Qn, aPolyhedron);
                        aRight->addEdge(Qm, Qn, aPolyhedron);
                        }
                        else
                        {
                            V_ASSERT(0);
                        }
                    }
                }
            }
        }

        aLeft->removeCollapsedEdges(aPolyhedron,tolerance);
        aRight->removeCollapsedEdges(aPolyhedron,tolerance);
        V_ASSERT(aLeft->getEdgeCount() > 0);
        V_ASSERT(aRight->getEdgeCount() > 0);
        //V_ASSERT(aLeft->isValid(aPolyhedron, normalization,tolerance));
        //V_ASSERT(aRight->isValid(aPolyhedron, normalization,tolerance));
        
        return ON_BOUNDARY;
    }

template<class S>
class SplitAlgorithmStatus
{
public:
    GeometryPositionType getPosition(PluckerElement* aElement) const
    {
        auto myIter = mPositions.find(aElement);
        if (myIter != mPositions.end())
        {
            return myIter->second;
        }
        return GeometryPositionType::ON_UNKNOWN_POSITION;
    }

    void setPosition(PluckerElement* aElement, GeometryPositionType aPosition)
    {
        mPositions[aElement] = aPosition;
    }

    void setUnchecked(PluckerElement* aElement)
    {
        mUnchecked[aElement] = true;
    }

    void setChecked(PluckerElement* aElement)
    {
        mUnchecked[aElement] = false;
    }

    bool isUnchecked(PluckerElement* aElement) const
    {
        auto myIter = mUnchecked.find(aElement);
        if (myIter != mUnchecked.end())
        {
            return myIter->second;
        }
        return false;
    }

    void setPlaneOffset(PluckerVertex<P>* v, S result)
    {
        mPlaneOffsets[v] = result;
    }
    
    S getPosition(PluckerVertex<P>* v) const
    {
        auto myIter = mPlaneOffsets.find(v);
        if (myIter != mPlaneOffsets.end())
        {
            return myIter->second;
        }
        return 0;
    }
    
    void clear()
    {
        mPositions.clear();
        mUnchecked.clear();
    }
private:
        std::map<PluckerElement*, GeometryPositionType> mPositions;
        std::map<PluckerElement*, bool> mUnchecked;
        std::map<PluckerVertex<P>*, S> mPlaneOffsets;
};

template<class P, class S>
inline void split(const P& aPlane, PluckerPolytopeComplex<P>& aComplex, SplitAlgorithmStatus<S>& aStatus)
{           
    for (auto iter = aComplex.beginVertices(); iter != aComplex.endVertices(); iter++)
    {
        PluckerVertex<P>* v = *iter;
        S result = aPlane.dot(v.getPlucker());
        GeometryPositionType position = MathPredicates::getRelativePosition(S, tolerance);
        aStatus.setPlaneOffset(v,result);
        aStatus.setPosition(v,position);
    }
    
    std::vector<PluckerInterpolatedVertex*> newVertices;
    
    for (int k = 1; k < aComplex.getDimension(); k++)
    {
        std::vector<PluckerElement*> newElements;
        std::vector<PluckerElement*> newChildren;
        
        for (auto c: aComplex.getElements(k))
        {
            GeometryPositionType myPositionC = c->getChildrenPosition();
            aStatus.setPosition(c, myPositionC); 
            if (myPositionC == ON_BOUNDARY)
            {
                PluckerElement* f = NULL;
                
                if (k==EntityDimension::EDGE) // edges
                {
                    PluckerEdge<P>* myEdge = static_cast<PluckerEdge<P>*>(c);
                    PluckerInterpolatedVertex* v = new PluckerInterpolatedVertex(myEdge, 
                                                                                 aStatus.getPosition(myEdge->getVertex0()), 
                                                                                 aStatus.getPosition(myEdge->getVertex1()));
                    newVertices.push_back(v);
                    f = v;
                 }
                 else
                 {
                    f = new PluckerElement();
                 }
                newChildren.push_back(f);
                aStatus.setPosition(f, ON_BOUNDARY);
               
                for (auto km2: c->getChildren())
                {
                   if (aStatus.getPosition(km2) == ON_BOUNDARY)
                   {
                       f->appendChildren(km2);
                   }
                }

                PluckerElement* cp = new PluckerElement();
                PluckerElement* cm = new PluckerElement();
                newChildren.push_back(cp);
                newChildren.push_back(cm);

                cp->appendChildren(f); 
                cm->appendChildren(f);
                
                for (auto kp1: c->getParents())
                {
                    cp->appendParent(f);
                    cm->appendParent(f);
                }

                for (auto km1: c->getChildren())
                {
                    GeometryPositionType myPositionKm1 = aStatus.getPosition(km1);
                    if (myPositionKm1== ON_NEGATIVE_SIDE)
                    {
                        cm->appendParent(f);
                    }                                              
                    else if (myPositionKm1 == ON_POSITIVE_SIDE)                    
                    {
                       cp->appendParent(f);
                    }
                }

                aComplex.deleteElement(c,k);

                for (auto element: newChildren)
                {
                    aComplex.appendElement(element,k-1);                    

                }
                for (auto element: newElements)
                {
                    aComplex.appendElement(element,k);
                }              
            }
        }

    }

    for (auto vertex: newVertices)
    {
        vertex.interpolate();
    }
}

template<class P, class S>
inline void reclassify(PluckerPolytopeComplex<P>& aComplex, SplitAlgorithmStatus<S>& aStatus)
{
    //Mark all the polytopes that have at least one vertex classified ON_BOUNDARY as unchecked
    
    for (auto iter = aComplex.beginVertices(); iter != aComplex.endVertices(); iter++)
    {
        PluckerVertex<P>* vertex = *iter;
    
        if (aStatus.getPosition(vertex) == ON_BOUNDARY)
        {
            std::stack<PluckerElement*> stack;
            for (auto parent: vertex->getParents())
            {
                stack.push(parent);
            }
            while (!stack.empty())
            {
                PluckerElement* element = stack.top();
                stack.pop();

                aStatus.setChecked(element);
                for (auto parent: element->getParents())
                {
                    stack.push(parent);
                }                
            }
        }        
    }

    for (int k = 1; k < myLattice->getDimension(); k++)
    {
        //for each polytope c of intrinsic dimension k ...
        for (auto c: myLattice->getElements(k))
        {
            // ... that is unchecked
            if (aStatus.isUnchecked(c))
            {
                // Mark c checked
                aStatus.setChecked(c);
                GeometryPositionType position = c->getChildrenPosition();
                aStatus.setPosition(c,position);

                if (false)
                {
                    //classify all the vertices of c as ON_BOUNDARY
                    std::stack<PluckerElement*> stack;
                    for (auto child: c->getChildren())
                    {
                        stack.push(child);
                    }
                    while (!stack.empty())
                    {
                        PluckerElement* element = stack.top();
                        stack.pop();  
                        if (element->getRank() == 0)
                        {
                            aStatus.setPosition(element, ON_BOUNDARY);
                        }
                        else
                        {
                            for (auto child: element->getChildren())
                            {
                                stack.push(child);
                            }  
                        }
                    }

                    // Mark all the faces of c as unchecked
                    for (auto child: c->getChildren())
                    {
                        aStatus.setUnchecked(child)
                    }

                    //set k=0 and continue to step 2 (we set k =1 for to re-enter the for loop)
                    k=1;
                    break;
                }
            }
        }
    }
}
}