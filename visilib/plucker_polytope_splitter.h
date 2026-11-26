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
#include "geometry_position_type.h"
#include "math_plucker_normalization.h"
#include "math_predicates.h"
#include "math_geometry.h"
#include "math_combinatorial.h"
#include "plucker_polytope.h"
#include "plucker_polyhedron.h"

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
        @param aNormalization: how the new vertices are normalized
        @param tolerance: a tolerance to determine  wheter a point lies on the splitting hyperplane or not
        @param aLeft: the resulting splitted polytope at the negative side of the hyperplane (must be instanciated before calling the function)
        @param aRight: the resulting splitted polytope at the positive side of the hyperplane  (must be instanciated before calling the function)
        @param aStatistic: where the incidences and the merges are counted (optional)
       */
        static GeometryPositionType split(PluckerPolyhedron<P>* polyhedron, const P& aPlane, PluckerPolytope<P>* aPolytope, PluckerPolytope<P>* aLeft, PluckerPolytope<P>* aRight, size_t aPlaneId, const MathPluckerNormalization<P, S>& aNormalization, S tolerance, HelperStatisticCollector* aStatistic = nullptr);

    private:
        /** @brief Whether two vertices are the same point: within the tolerance, or, with the symbolic perturbation, equal
        in every term of epsilon (an edge only epsilon long joins two different points and is kept).*/
        static bool areSameVertex(const PluckerPolyhedron<P>* aPolyhedron, size_t a, size_t b, S tolerance)
        {
            if constexpr (MathSoSSupported<P, S>::value)
            {
                if (aPolyhedron->isSymbolic())
                {
                    return MathSoS<P, S>::areEqual(aPolyhedron->getSymbolicPoint(a), aPolyhedron->getSymbolicPoint(b));
                }
            }
            return MathPredicates::isEdgeCollapsed(aPolyhedron->get(a), aPolyhedron->get(b), tolerance);
        }

        /** @brief Whether the polyhedron carries the symbolic perturbation (only ever for an exact arithmetic)*/
        static bool isSymbolic(const PluckerPolyhedron<P>* aPolyhedron)
        {
            if constexpr (MathSoSSupported<P, S>::value)
            {
                return aPolyhedron->isSymbolic();
            }
            return false;
        }

        /** @brief The side of the vertex aVertex of the perturbed hyperplane aPlaneID (see MathSoS::sign)*/
        static GeometryPositionType symbolicSide(const PluckerPolyhedron<P>* aPolyhedron, const P& aPlane, size_t aPlaneID, size_t aVertex, HelperStatisticCollector* aStatistic)
        {
            if constexpr (MathSoSSupported<P, S>::value)
            {
                // the order 0 term is the stored point: it decides almost always, without building the symbolic point
                const int leading = getExactSign(aPlane.dot(aPolyhedron->get(aVertex)));
                if (leading != 0)
                {
                    return leading > 0 ? ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE;
                }
                size_t order = 0;
                const int s = MathSoS<P, S>::sign(aPlane, aPolyhedron->getPerturbation(aPlaneID), aPolyhedron->getSymbolicPoint(aVertex), order);
                if (aStatistic != nullptr)
                {
                    if (s == 0)
                    {
                        aStatistic->inc(SOS_RESIDUAL_COUNT);
                    }
                    else if (order > 0)
                    {
                        aStatistic->inc(SOS_TIE_BREAK_COUNT);
                    }
                }
                return s > 0 ? ON_POSITIVE_SIDE : (s < 0 ? ON_NEGATIVE_SIDE : ON_BOUNDARY);
            }
            return ON_BOUNDARY;
        }

        /** @brief Build, normalize and add the vertex where the edge from aNegative to aPositive crosses the perturbed
        hyperplane aPlaneID; its side of the quadric is the sign of its symbolic quadric. Return its index.*/
        static size_t addSymbolicCrossing(PluckerPolyhedron<P>* aPolyhedron, const P& aPlane, size_t aPlaneID, size_t aNegative, size_t aPositive,
                                          const MathPluckerNormalization<P, S>& aNormalization, S tolerance)
        {
            if constexpr (MathSoSSupported<P, S>::value)
            {
                MathSymbolicPoint<P> crossing = MathSoS<P, S>::crossing(aPlane, aPolyhedron->getPerturbation(aPlaneID),
                                                                         aPolyhedron->getSymbolicPoint(aNegative), aPolyhedron->getSymbolicPoint(aPositive));
                aNormalization.normalizeSymbolic(crossing);
                size_t order = 0;
                const int q = MathSoS<P, S>::quadricSign(crossing, order);
                const GeometryPositionType quadric = q > 0 ? ON_POSITIVE_SIDE : (q < 0 ? ON_NEGATIVE_SIDE : ON_BOUNDARY);
                const size_t index = aPolyhedron->add(crossing.leading(), quadric, aNormalization.isUnitNorm(), tolerance);
                aPolyhedron->setSymbolicPoint(index, crossing);
                return index;
            }
            return 0;
        }
    };

    template<class P, class S>
    inline GeometryPositionType PluckerPolytopeSplitter<P, S>::split(PluckerPolyhedron<P>* aPolyhedron, const P& aPlane, PluckerPolytope<P>* aPolytope, PluckerPolytope<P>* aLeft, PluckerPolytope<P>* aRight, size_t aPlaneID, const MathPluckerNormalization<P, S>& normalization, S tolerance, HelperStatisticCollector* aStatistic)
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


        // Step 1 - edge classification

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

                    /* A vertex the splitting hyperplane is already a facet of lies on it by construction: the answer is
                    read from its facet description, and the rounded value of the product, which may stand a little
                    off zero, is never used for it. Only the other vertices are classified numerically. */
                    S result(0);
                    GeometryPositionType position = ON_BOUNDARY;
                    if (MathCombinatorial::hasFacet(aPolyhedron->getFacetsDescription(a), aPlaneID))
                    {
                        if (aStatistic != nullptr)
                        {
                            aStatistic->inc(INCIDENCE_COMBINATORIAL_COUNT);
                        }
                    }
                    else if (isSymbolic(aPolyhedron))
                    {
                        /* Simulation of Simplicity: the side is the sign of the polynomial of epsilon of the perturbed
                        hyperplane at the symbolic vertex. A polynomial vanishing to every order kept is an incidence by
                        construction, answered like the recorded ones. */
                        position = symbolicSide(aPolyhedron, aPlane, aPlaneID, a, aStatistic);
                    }
                    else
                    {
                        result = aPlane.dot(aPolyhedron->get(a));

                        /* A vertex and a hyperplane that are both lines through grid points: the side is the exact sign of
                        their Klein product, and a zero is a true incidence. The rounded value is kept for the
                        interpolation only. The other vertices are classified numerically. */
                        int exactSign = 0;
                        if (aPolyhedron->getExactSign(aPlaneID, a, exactSign))
                        {
                            position = exactSign > 0 ? ON_POSITIVE_SIDE : (exactSign < 0 ? ON_NEGATIVE_SIDE : ON_BOUNDARY);
                            if (aStatistic != nullptr)
                            {
                                aStatistic->inc(INCIDENCE_EXACT_COUNT);
                            }
                        }
                        else
                        {
                            position = MathPredicates::getRelativePosition(result, tolerance);
                            V_ASSERT(position == MathPredicates::getVertexPlaneRelativePosition(aPlane, aPolyhedron->get(a), tolerance));
                            if (position == ON_BOUNDARY && aStatistic != nullptr)
                            {
                                aStatistic->inc(INCIDENCE_NUMERIC_BOUNDARY_COUNT);
                            }
                        }
                    }

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
                                                  myFacets,
                                                  MathPluckerCellTraits<P>::cCommonFacetsForEdge);


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
                    size_t vertexIndex = 0;
                    if (isSymbolic(aPolyhedron))
                    {
                        // the exact crossing of the perturbed hyperplane, carrying its terms of epsilon (myI1 is the negative end)
                        vertexIndex = addSymbolicCrossing(aPolyhedron, aPlane, aPlaneID, myI1, myI2, normalization, tolerance);
                    }
                    else
                    {
                        P  myIntersection = MathGeometry::interpolate(myArrayOfPlanePosition[myI1], myArrayOfPlanePosition[myI2], aPolyhedron->get(myI1), aPolyhedron->get(myI2), tolerance);
                        myIntersection = normalization.normalizePoint(myIntersection);

                        vertexIndex = aPolyhedron->add(myIntersection, MathPredicates::getQuadricRelativePosition(myIntersection, tolerance), normalization.isUnitNorm(), tolerance);
                    }

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

                if (MathCombinatorial::haveAtLeastNCommonFacets(facetsQm, facetsQn, MathPluckerCellTraits<P>::cCommonFacetsForEdge))
                {
                    V_ASSERT(Qm != Qn);
                    if (!areSameVertex(aPolyhedron, Qn, Qm, tolerance))
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

        aLeft->removeCollapsedEdges(aPolyhedron, tolerance, aStatistic);
        aRight->removeCollapsedEdges(aPolyhedron, tolerance, aStatistic);
        V_ASSERT(aLeft->getEdgeCount() > 0);
        V_ASSERT(aRight->getEdgeCount() > 0);
        //V_ASSERT(aLeft->isValid(aPolyhedron, normalization.isUnitNorm(), tolerance));
        //V_ASSERT(aRight->isValid(aPolyhedron, normalization.isUnitNorm(), tolerance));

        return ON_BOTH_SIDES;
    }
}
