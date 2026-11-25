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

#include <set>
#include <unordered_set>
#include "math_plucker_6.h"
#include "math_geometry.h"
#include "math_combinatorial.h"
#include "geometry_position_type.h"
#include "plucker_polyhedron.h"

namespace visilib
{
    template<class P>
    class PluckerPolyhedron;

    class Silhouette;

    /** @brief Store a polytope in Plucker space.

    Only the skeleton (vertices and edges) of the polytope is encoded .*/

    template<class P>
    class PluckerPolytope
    {
    public:
        PluckerPolytope<P>();

        ~PluckerPolytope<P>();

        /** @brief Return an extremal stabbing line ESL of the polytope.

        The ESL lies at the intersection of the edges of the polytope with the Plucker quadric
        @param i: index of the ESL
        @return: the ESL equation
        */
        inline const P getExtremalStabbingLine(size_t i)
        {
            return mExtremalStabbingLines[i];
        }


        /** @brief Return the number of extremal stabbing lines */
        inline size_t getExtremalStabbingLinesCount() const
        {
            return mExtremalStabbingLines.size();
        }


        /** @brief Return the number of edges of the polytope */

        inline size_t getEdgeCount() const
        {
            return mEdges.size();
        }

        /** @brief Return the edges of the polytope, wich are hyperlines joining two vertices
        @return: a set containing all the edges <i,j> joining the vertices i and j
        */
        inline const std::set<std::pair<size_t, size_t>>& getEdges() const
        {
            return mEdges;
        }


        /** @brief Return the number of vertices of the polytope */
        inline const std::unordered_set<size_t>& getVertices() const
        {
            return mVertices;
        }

        /** @brief Return the subset of edges of the polytope that have an intersection with the Plucker Quadric (ie reprsents real line)
       @return: a set containing all the edges <i,j> joining the vertices i and j with an intersection with the Plucker Quadric
       */

        inline const std::set<std::pair<size_t, size_t>>& getEdgesIntersectingQuadric()
        {
            return mEdgesIntersectingQuadric;
        }

        /** @brief Add an edge to the polytope
        @param aVertex0: the first vertex of the edge
        @param aVertex1: the second vertex of the edge
        @param aPolyhedron: the polyhedron that contains the vertices equation
        */
        void addEdge(size_t aVertex0, size_t aVertex1, PluckerPolyhedron<P>* aPolyhedron);


        /** @brief Test if the polytope is valid (fast check)

        A polytope is containsRealLines if it has edges, if the number of edges is not too high and if at least one edge crosses the Plucker Quadric
        */
        inline bool containsRealLines() const;

        /** @brief Test if the polytope is degenerated (the length of all edges are below the given tolerance).

        @param polyhedron: the polyhedron containing all the vertices definition
        @param S: the tolerance against which the length of the edges is tested
        */

        template<class S> bool hasSomeEdgesCollapsed(PluckerPolyhedron<P>* polyhedron, S Tolerance);

        /** @brief Remove the collasped edge (the length of all edges are below the given tolerance).

        @param polyhedron: the polyhedron containing all the vertices definition
        @param S: the tolerance against which the length of the edges is tested
        */

        template<class S> bool removeCollapsedEdges(PluckerPolyhedron<P>* polyhedron, S Tolerance, HelperStatisticCollector* aStatistic = nullptr);


        /** @brief Compute the extremal stabbing lines of the polytope*/
        template<class S> void computeExtremalStabbingLines(PluckerPolyhedron<P>* polyhedron, S Tolerance);


        /** @brief Compute the list of edges intersecting the Plucker quadric*/
        template<class S> void computeEdgesIntersectingQuadric(PluckerPolyhedron<P>* aPolyhedron, S Tolerance);


        /** @brief Convert all the extremal stabbing lines to 3D space.

        The intersection of each extremal stabbing line with two 3D planes is computed, the two resulting points define the line in 3D.
        @param aStabbingLines: vector containing pair of points defining the stabbing lines in 3D. Each line is represented by two points
        @param aPlane1: first plane that will be intersected by each stabbing lines to define the first point of the extremal line in 3D
        @param aPlane2: second plane that will be intersected by each stabbing lines to define the second point of the extremal line in 3D
      */
        void getExtremalStabbingLinesBackTo3D(std::vector<std::pair<MathVector3d, MathVector3d> >& aStabbingLines, const MathPlane3d& aPlane1, const MathPlane3d& aPlane2);


        /** @brief Add an extremal stabbing line to the polytope.

        The ESL lies at the intersection of the edges of the polytope with the Plucker quadric
        */
        void addExtremalStabbingLine(const P& aStabbingLine)
        {
            mExtremalStabbingLines.push_back(aStabbingLine);
        }


        /** @brief Output in the console the properties of the polytope.*/
        void outputProperties(std::ostream& out, PluckerPolyhedron<P>* polyhedron = nullptr);


        /** @brief Returns the list of all the faces (ie the equation of the hyperplane) defining the polytope.

        The hypeplanes are stored in a polyhedon , and their facets is represented as an index list inside the polyhedron
        */
        void getFacets(std::set<size_t>& facets, PluckerPolyhedron<P>* polyhedron);

        /** @brief  One facet of the polytope, in double precision, with the side the polytope lies on

        The plane is divided by its norm, so that the value it takes on a point of unit norm is a distance.*/
        struct CellFacet
        {
            MathVector3d direction;   /**< @brief The direction part of the plane*/
            MathVector3d location;    /**< @brief The location part of the plane*/
            double sign;              /**< @brief 1 when the cell lies where the plane is negative, -1 where it is positive, 0 when the whole cell lies in the plane*/
            size_t index;             /**< @brief The index of the plane in the polyhedron*/
        };

        /** @brief  Collect the facets of the polytope, that is the planes tight on one of its vertices

        Every facet of the convex hull of the vertices carries one of these planes, so they describe the cell
        completely. A plane the whole cell lies in comes back with a sign of nought, an equality a line of the cell
        has to lie in as well, and is what the flat flag tells of: behind such a plane there is no interior to reach.*/
        template<class S> size_t getCellFacets(PluckerPolyhedron<P>* polyhedron, const S & tolerance, std::vector<CellFacet>& facets, bool& flat);

        /** @brief  Collect, among the planes of someCandidates, those the polytope lies on one side of (or in), with that side

        The planes tight on a vertex are not all the planes the cell lies behind: a plane that cut an ancestor of the cell,
        or an edge of a source, can leave every vertex strictly on one side, and the hyperplanes the cell records are only a
        part of the facets of the convex hull of its vertices in the six dimensional space. Their intersection holds lines
        the hull does not; any plane the cell lies behind is a valid bound of the hull, and adding them gives back the real
        lines of the cell. The planes with vertices on both sides are left out.*/
        template<class S, class I> size_t getCellFacets(PluckerPolyhedron<P>* polyhedron, const S & tolerance, const I& someCandidates, std::vector<CellFacet>& facets, bool& flat);

        /** @brief Add a silhouette to the polyptope */
        void addSilhouette(Silhouette* silhouette)
        {
            if (mSilhouettes.find(silhouette) == mSilhouettes.end())
            {
                mSilhouettes.insert(silhouette);
            }
        }

        void setRepresentativeLine(P line)
        {
            mRepresentativeLine = line;
        }

        P getRepresentativeLine() const
        {
            return mRepresentativeLine;
        }

        void setRepresentativeLineMargin(double margin)
        {
            mRepresentativeLineMargin = margin;
        }

        double getRepresentativeLineMargin() const
        {
            return mRepresentativeLineMargin;
        }

        /** @brief Test if the polytope is valid (complete check)

        The test performs a complete validation of the validity of the polytope
        */
        template<class S> bool isValid(PluckerPolyhedron<P>* polyhedron, bool normalization, S tolerance);

    private:

        std::set<std::pair<size_t, size_t> > mEdges;                      /** < @brief The edges of the polytope, each edge is defined by the indices of the two vertices supporting it*/
        std::set<std::pair<size_t, size_t> > mEdgesProcessed;
        std::vector<P> mExtremalStabbingLines;                           /** < @brief The ESL of the polytope, at the intersection of an edge and the Plucker Quadric*/
        std::set<std::pair<size_t, size_t> > mEdgesIntersectingQuadric;   /** < @brief The edges containing an intersection with the Plucker Quadric.*/
        std::unordered_set<Silhouette*> mSilhouettes;          /** < @brief The set of silhouettes associated to the polytope*/
        std::unordered_set<size_t> mVertices;                            /** < @brief The indices of the vertices of the polytope*/
        std::vector<std::vector<size_t> > mExtremalStabbingLinesFacets;
        double mRadius;
        P mRepresentativeLine;
        /** @brief How far the line above stands from the nearest facet of the polytope, in its own arithmetic

        A representative line is what the search casts when it wants to know about a cell, and the margin is
        what says how much room the cast has: a line of no margin was accepted as good as by luck, and is the
        one to look at when the cast and the cell disagree. */
        double mRepresentativeLineMargin;
    };

    template<class P>
    PluckerPolytope<P>::PluckerPolytope()
        :    mRadius(0),
             mRepresentativeLineMargin(0.0)
    {
    }

    template<class P>
    PluckerPolytope<P>::~PluckerPolytope()
    {
    }

    template<class P>
    void PluckerPolytope<P>::addEdge(size_t aVertex0, size_t aVertex1, PluckerPolyhedron<P>* aPolyhedron)
    {
        size_t min = aVertex0;
        size_t max = aVertex1;
        if (aVertex0 > aVertex1)
        {
            min = aVertex1;
            max = aVertex0;
        }
        V_ASSERT(aVertex0 != aVertex1);
        V_ASSERT(MathCombinatorial::haveAtLeastNCommonFacets(aPolyhedron->getFacetsDescription(aVertex0), aPolyhedron->getFacetsDescription(aVertex1), MathPluckerCellTraits<P>::cCommonFacetsForEdge));
        if (mEdges.find(std::pair<size_t, size_t>(min, max)) == mEdges.end())
        {
            mEdges.insert(std::pair<size_t, size_t>(min, max));
        }

        if (mVertices.find(aVertex0) == mVertices.end())
        {
            mVertices.insert(aVertex0);
        }
        if (mVertices.find(aVertex1) == mVertices.end())
        {
            mVertices.insert(aVertex1);
        }
    }

    template<class P>
    void PluckerPolytope<P>::getExtremalStabbingLinesBackTo3D(std::vector<std::pair<MathVector3d, MathVector3d>>& aStabbingLines, const MathPlane3d& aPlane1, const MathPlane3d& aPlane2)
    {
        aStabbingLines.clear();
        V_ASSERT(mExtremalStabbingLines.size() > 0);
        for (size_t i = 0; i < mExtremalStabbingLines.size(); i++)
        {
            std::pair<MathVector3d, MathVector3d> r = MathGeometry::getBackTo3D(mExtremalStabbingLines[i], aPlane1, aPlane2);
            aStabbingLines.push_back(r);
        }

#if 1

        if (mExtremalStabbingLinesFacets.size() == 0)
            return;
        for (size_t i = 0; i < mExtremalStabbingLines.size(); i++)
        {
            for (size_t j = i + 1; j < mExtremalStabbingLines.size(); j++)
            {
                if (MathCombinatorial::haveAtLeastNCommonFacets(mExtremalStabbingLinesFacets[i], mExtremalStabbingLinesFacets[j], MathPluckerCellTraits<P>::cCommonFacetsForEdge))
                {
                    aStabbingLines.push_back(std::pair<MathVector3d, MathVector3d>( aStabbingLines[i].first,  aStabbingLines[j].first));
                    aStabbingLines.push_back(std::pair<MathVector3d, MathVector3d>( aStabbingLines[i].second, aStabbingLines[j].second));
                }
            }
        }
#endif

        V_ASSERT(aStabbingLines.size() > 0);
    }

    template<class P>
    void PluckerPolytope<P>::outputProperties(std::ostream& o, PluckerPolyhedron<P>* polyhedron)
    {
        o << "Polytope ESL: " << mExtremalStabbingLines.size() << std::endl;
        o << "Polytope Edges: " << mEdges.size() << std::endl;
        o << "Polytope Vertices: " << mVertices.size() << std::endl;
        o << "Polytope representative line margin: " << mRepresentativeLineMargin << std::endl;

        if (polyhedron != nullptr)
        {
            o << "Polyhedron Lines:" << polyhedron->getLinesCount() << std::endl;
            std::set<size_t> myFacets;
            getFacets(myFacets, polyhedron);

            o << "Facets: {";

            for (auto i = myFacets.begin(); i != myFacets.end(); i++)
            {
                o << *i << ", ";
    }
            o << "}" << std::endl;

            o << "Polytope vertices facets: " << std::endl;
            for (auto iter = mVertices.begin(); iter != mVertices.end(); iter++)
            {
                const std::vector<size_t>& facets = polyhedron->getFacetsDescription(*iter);
                o << " v[" << *iter << "] : facets{";
                for (size_t i = 0; i < facets.size(); i++)
                {
                    o << facets[i] << ", ";
                }
                o << "}";

                o << polyhedron->get(*iter);
                o << std::endl;
            }

        }
            o << "Polytope edges facets: " << std::endl;

            int line = 0;
            for (auto iter = mEdges.begin(); iter != mEdges.end(); iter++)
            {
                o << " e[" << iter->first << "," << iter->second << "]";
            }
            o << std::endl;

        }
    template<class P> template<class S>
    bool PluckerPolytope<P>::removeCollapsedEdges(PluckerPolyhedron<P>* polyhedron, S tolerance, HelperStatisticCollector* aStatistic)
    {
        std::vector<int> myEdgesTable;
        std::vector<int> myMergeTable;

        for (auto iter = mEdges.begin(); iter != mEdges.end();iter++)
        {
            size_t i1 = iter->first;
            size_t i2 = iter->second;

            const P& v1 = polyhedron->get(i1);
            const P& v2 = polyhedron->get(i2);


            if (MathPredicates::isEdgeCollapsed(v1, v2, tolerance))
            {
                V_TRACE("collapsed edge dected: " << v1-v2 << "; tolerance: "<< 0.5 *tolerance);
                myMergeTable.push_back(i1);
                myMergeTable.push_back(i2);
            }
            else
            {
                myEdgesTable.push_back(i1);
                myEdgesTable.push_back(i2);
            }
        }

        if (myMergeTable.empty())
        {
            return false;
        }

#ifdef OUTPUT_DEBUG_FILE
        std::cout << "EDGE TO BE MERGED DETECTED!" << std::endl;
        outputProperties(std::cout,polyhedron);

        std::cout << "myMergeTable table:";
        for (int i=0; i<myMergeTable.size();i++)
        {
            std::cout<<myMergeTable[i]<<" ";
        }
        std::cout<<std::endl;
        std::cout << "myEdgesTable table:";
        for (int i=0; i<myEdgesTable.size();i++)
        {
            std::cout<<myEdgesTable[i]<<" ";
        }
        std::cout<<std::endl;
#endif
        for (int i=0; i<myMergeTable.size();i+=2)
        {
            int merge_left  = myMergeTable[i];
            int merge_right = myMergeTable[i+1];
            if (merge_left!=merge_right)
            {
                if (aStatistic != nullptr)
                {
                    aStatistic->inc(COLLAPSED_EDGE_MERGED_COUNT);
                }
                for (int j=0; j<myEdgesTable.size();j++)
                {
                    if (myEdgesTable[j] == merge_right)
                    {
                        myEdgesTable[j] = merge_left;
                    }
                }

                for (int k=i;k<myMergeTable.size();k++)
                {
                    if (myMergeTable[k] == merge_right)
                    {
                        myMergeTable[k] = merge_left;
                    }
                }

                const std::vector<size_t>& facetsI1 = polyhedron->getFacetsDescription(merge_left);
                const std::vector<size_t>& facetsI2 = polyhedron->getFacetsDescription(merge_right);

                std::vector<size_t> myFacets;
                MathCombinatorial::initFacets(facetsI1,
                                              facetsI2,
                                              myFacets);
                polyhedron->initFacetsDescription(merge_left, myFacets);
            }
        }

#ifdef OUTPUT_DEBUG_FILE
        std::cout << "myMergeTable table:";
        for (int i=0; i<myMergeTable.size();i++)
        {
            std::cout<<myMergeTable[i]<<" ";
        }
        std::cout<<std::endl;
        std::cout << "myEdgesTable table:";
        for (int i=0; i<myEdgesTable.size();i++)
        {
            std::cout<<myEdgesTable[i]<<" ";
        }
        std::cout<<std::endl;

        outputProperties(std::cout,polyhedron);
#endif
        mEdges.clear();
        for (int i=0; i<myEdgesTable.size();i+=2)
        {
            addEdge(myEdgesTable[i], myEdgesTable[i+1],polyhedron);
        }
	return true;
    }

    template<class P> template<class S>
    bool PluckerPolytope<P>::hasSomeEdgesCollapsed(PluckerPolyhedron<P>* polyhedron, S tolerance)
    {
        for (auto iter = mEdges.begin(); iter != mEdges.end(); iter++)
        {
            const P& v1 = polyhedron->get(iter->first);
            const P& v2 = polyhedron->get(iter->second);

            if (MathPredicates::isEdgeCollapsed(v1, v2, tolerance))
            {
                return true;
            }
            const std::vector<size_t>& facets1 = polyhedron->getFacetsDescription(iter->first);
            const std::vector<size_t>& facets2 = polyhedron->getFacetsDescription(iter->second);

            if (facets1.size() == facets2.size() && std::equal(facets1.begin(), facets1.end(), facets2.begin()))
            {
                return true;
            }
        }

        return false;
    }

    template<class P> template<class S>
    void PluckerPolytope<P>::computeExtremalStabbingLines(PluckerPolyhedron<P>* polyhedron, S tolerance)
    {
        V_ASSERT(mExtremalStabbingLines.empty());
        V_ASSERT(!mEdgesIntersectingQuadric.empty());

        for (auto iter = mEdgesIntersectingQuadric.begin(); iter != mEdgesIntersectingQuadric.end(); iter++)
        {
            const P& v1 = polyhedron->get(iter->first);
            const P& v2 = polyhedron->get(iter->second);
            GeometryPositionType p1 = polyhedron->getQuadricRelativePosition(iter->first);
            GeometryPositionType p2 = polyhedron->getQuadricRelativePosition(iter->second);

            std::vector<P> p;

            const std::vector<size_t>& facets1 = polyhedron->getFacetsDescription(iter->first);
            const std::vector<size_t>& facets2 = polyhedron->getFacetsDescription(iter->second);
            std::vector<size_t> edgeFacets;
            MathCombinatorial::initFacets(facets1, facets2, polyhedron->getLinesCount(), edgeFacets, MathPluckerCellTraits<P>::cCommonFacetsForEdge);
            edgeFacets.pop_back();

            if (MathGeometry::findPluckerEdgeWithQuadricIntersection<P, S>(v1, v2, p1, p2, p, tolerance))
            {
                for (size_t index = 0; index < p.size(); index++)
                {
                    //	V_ASSERT(MathPredicates<P>::isNormalized(p[index]));
                    addExtremalStabbingLine(p[index]);

                    mExtremalStabbingLinesFacets.push_back(edgeFacets);
                }
            }
            else
            {
                V_ASSERT(MathGeometry::hasPluckerEdgeWithQuadricIntersection(v1, v2, p1, p2, tolerance));
            }
        }
    }

    template<class P>
    bool PluckerPolytope<P>::containsRealLines() const
    {
        if (mEdgesIntersectingQuadric.size() == 0)
        {
            return false;
        }

        if (mEdges.size() == 0)
        {
            return false;
        }
        if (mEdges.size() > 5000)
        {
            V_ASSERT(0);
            return false;
        }
        V_ASSERT(mVertices.size() > 0);
        return true;
    }

    template<class P> template<class S>
    bool PluckerPolytope<P>::isValid(PluckerPolyhedron<P>* polyhedron, bool normalization, S tolerance)
    {
        bool isValid = true;

        if (mEdgesIntersectingQuadric.size() == 0)
        {
            std::cerr << "   IsValid error: no stabbing line..." << std::endl;
            isValid = false;
        }

        if (mEdges.size() == 0)
        {
            std::cerr << "   IsValid error: no edge..." << std::endl;

            isValid = false;
        }

        size_t edgeMax = 2000;
        if ((size_t)mEdges.size() > edgeMax)
        {
            std::cerr << "   IsValid error: more than " << edgeMax << " edges..." << std::endl;

            isValid = false;
        }

        if (MathPredicates::hasPluckerPolytopeIntersectionWithQuadric(this, polyhedron) != ON_BOUNDARY)
        {
            std::cerr << "   IsValid error: polytope not intersecting the quadric" << std::endl;

            isValid = false;
        }

        if (hasSomeEdgesCollapsed(polyhedron, tolerance))
        {
            std::cerr << "   IsValid error: hasSomeEdgesCollapsed" << std::endl;
            outputProperties(std::cout);

            isValid = false;
        }

        for (auto iter = mVertices.begin(); iter != mVertices.end(); iter++)
        {
            if (normalization)
            {
                if (!MathPredicates::isNormalized(polyhedron->get(*iter), tolerance))
                {
                    std::cerr << "   IsValid error: vertex " << *iter << " not normalized" << std::endl;

                    isValid = false;
                }
            }
            std::vector<size_t>& myFacets = polyhedron->getFacetsDescription(*iter);

            for (auto i = myFacets.begin(); i != myFacets.end(); i++)
            {
                const P& v1 = polyhedron->get(*iter);
                const P& v2 = polyhedron->get(*i);

                if (MathPredicates::getVertexPlaneRelativePosition(v1, v2, tolerance) != ON_BOUNDARY)
                {
                    std::cerr << "   IsValid error: polytope vertex " << *iter << " is not on the facet " << *i << std::endl;
                    std::cerr << "        P:" << polyhedron->get(*iter) << std::endl;
                    std::cerr << "        Facet:" << polyhedron->get(*i) << std::endl;
                    std::cerr << "        Facet.dot(P):" <<  polyhedron->get(*iter).dot(polyhedron->get(*i)) << std::endl;
                    isValid = false;
                }
            }
        }

        if (!isValid)
        {
            outputProperties(std::cout,polyhedron);
        }

        return isValid;
    }

    template<class P> template<class S>
    void PluckerPolytope<P>::computeEdgesIntersectingQuadric(PluckerPolyhedron<P>* polyhedron, S tolerance)
    {
        //Only compute if the function has not been called already
        if (mEdgesIntersectingQuadric.empty())
        {
            //At least one
            if (mEdges.size() > 0 && MathPredicates::hasPluckerPolytopeIntersectionWithQuadric(this, polyhedron))
            {
                V_ASSERT(mEdgesIntersectingQuadric.empty());
                V_ASSERT(!mVertices.empty());

                for (auto iter = mEdges.begin(); iter != mEdges.end(); iter++)
                {
                    V_ASSERT(std::find(mVertices.begin(), mVertices.end(), iter->first) != mVertices.end());
                    V_ASSERT(std::find(mVertices.begin(), mVertices.end(), iter->second) != mVertices.end());
                
                    const P& v1 = polyhedron->get(iter->first);
                    const P& v2 = polyhedron->get(iter->second);
                    GeometryPositionType p1 = polyhedron->getQuadricRelativePosition(iter->first);
                    GeometryPositionType p2 = polyhedron->getQuadricRelativePosition(iter->second);
                    if (MathGeometry::hasPluckerEdgeWithQuadricIntersection(v1, v2, p1, p2, tolerance))
                    {
                        mEdgesIntersectingQuadric.insert(*iter);
                    }
                }
            }
            else
            {
                return;
            }
            if (!containsRealLines())
            {
                    mEdgesIntersectingQuadric.clear();         
                   return;
            }
            V_ASSERT(mEdgesIntersectingQuadric.size() > 0);
        }
    }
    template<class P>
    void PluckerPolytope<P>::getFacets(std::set<size_t>& facets, PluckerPolyhedron<P>* polyhedron)
    {
        for (auto iter = mVertices.begin(); iter != mVertices.end(); iter++)
        {
            const std::vector<size_t>& f = polyhedron->getFacetsDescription(*iter);
            facets.insert(f.begin(), f.end());
        }
    }

    template<class P> template<class S>
    size_t PluckerPolytope<P>::getCellFacets(PluckerPolyhedron<P>* polyhedron, const S & tolerance, std::vector<CellFacet>& facets, bool& flat)
    {
        std::set<size_t> myFacetIndices;
        getFacets(myFacetIndices, polyhedron);

        return getCellFacets(polyhedron, tolerance, myFacetIndices, facets, flat);
    }

    template<class P> template<class S, class I>
    size_t PluckerPolytope<P>::getCellFacets(PluckerPolyhedron<P>* polyhedron, const S & tolerance, const I& someCandidates, std::vector<CellFacet>& facets, bool& flat)
    {
        facets.clear();
        flat = false;

        for (auto index = someCandidates.begin(); index != someCandidates.end(); index++)
        {
            const P plane = polyhedron->get(*index);

            const GeometryPositionType position = MathPredicates::getRelativePosition(this, polyhedron, plane, tolerance, size_t(*index));

            /* The cell lies on one side of its own facets, and the side is read from its vertices. A plane the whole
            cell lies in is kept too, as an equality: a line that leaves it is outside the cell, however well it
            satisfies the rest. A plane with vertices on both sides is not a facet of the hull and is left out. */
            if (position == ON_BOTH_SIDES) continue;

            /* The plane read in its own coordinates: the conversion to MathPlucker6d
            carried the same parts through, and its scalar product with itself is the
            six dimensional norm they make, whichever space holds them. */
            const MathVector3d planeDirection = convert<MathVector3d>(plane.getDirection());
            const MathVector3d planeLocation = convert<MathVector3d>(plane.getLocation());
            const double norm = MathArithmetic<double>::getSqrt(planeDirection.dot(planeDirection) + planeLocation.dot(planeLocation));

            if (!(norm > 0.0)) continue;

            CellFacet facet;
            facet.index = *index;

            if (position == ON_BOUNDARY)
            {
                facet.sign = 0.0;
                flat = true;
            }
            else
            {
                facet.sign = position == ON_NEGATIVE_SIDE ? 1.0 : -1.0;
            }
            facet.direction = MathVector3d(planeDirection.x / norm, planeDirection.y / norm, planeDirection.z / norm);
            facet.location = MathVector3d(planeLocation.x / norm, planeLocation.y / norm, planeLocation.z / norm);

            facets.push_back(facet);
        }

        return facets.size();
    }
}
