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
    template<class P, class S>
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

                S getPlaneOffset(PluckerVertex<P>* v) const
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
                static void split(const P& aPlane, PluckerPolytopeComplex<P,S>& aComplex, SplitAlgorithmStatus<P,S>& aStatus, S aTolerance);
        };

    template<class P, class S>
        inline void PluckerPolytopeSplitter<P,S>::split(const P& aPlane, PluckerPolytopeComplex<P,S>& aComplex, SplitAlgorithmStatus<P,S>& aStatus, S aTolerance)
        { 
            for (auto iter = aComplex.beginVertices(); iter != aComplex.endVertices(); iter++)
            {
                PluckerVertex<P>* v = *iter;
                S result = aPlane.dot(v.getPlucker());
                GeometryPositionType position = MathPredicates::getRelativePosition(result, aTolerance);
                aStatus.setPlaneOffset(v,result);
                aStatus.setPosition(v,position);
            }

            std::vector<PluckerInterpolatedVertex<P,S>*> newVertices;

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
                        PluckerElement* f = nullptr;
                        if (k==EDGE) // edges
                        {
                            PluckerEdge<P>* myEdge = static_cast<PluckerEdge<P>*>(c);
                            PluckerInterpolatedVertex<P,S>* v = new PluckerInterpolatedVertex<P,S>(myEdge, 
                                    aStatus.getPosition(myEdge->getVertex0()), 
                                    aStatus.getPosition(myEdge->getVertex1()));
                            newVertices.push_back(v);
                            f = v;
                        }
                        else
                        {
                            f = PluckerElementFactory::createElement<P,S>(k-1);
;
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

                        PluckerElement* cp = PluckerElementFactory::createElement<P,S>(k-1);
                        PluckerElement* cm = PluckerElementFactory::createElement<P,S>(k-1);
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
        inline void PluckerPolytopeSplitter<P,S>::reclassify(PluckerPolytopeComplex<P,S>& aComplex, SplitAlgorithmStatus<S>& aStatus)
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

