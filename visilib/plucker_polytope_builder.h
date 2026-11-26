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
#include "math_plucker_normalization.h"
#include "plucker_polyhedron.h"
#include "visibility_event_registry.h"

namespace visilib
{
    /** @brief Create the minimal polytope containing all the stabbing lines for polygon a and b */

    template <class P, class S>
    class PluckerPolytopeBuilder
    {
    public:
        /** @brief A builder of the minimal polytope. When aRegistry is given, the edges of the two polygons are recorded in it as the
        features their hyperplanes carry: the first polygon as the light, the second as the receiver.*/
        PluckerPolytopeBuilder(const MathPluckerNormalization<P, S>& aNormalization, S tolerance, VisibilityEventRegistry* aRegistry = nullptr, const MathIntegerGrid& aGrid = MathIntegerGrid());

        /** @brief Create the minimal polytope following Mora et al. An extermal stabbing line is created for each vertex pair of polygon a and b.
        @param a: first source polygon
        @param b: second source polygon
        @param polyhedron: the polyhedron that will contains all the vertices of the polytope
        @return: the created polytope
        */
        PluckerPolytope<P>* build(GeometryConvexPolygon& a, GeometryConvexPolygon& b, PluckerPolyhedron<P>* polyhedron);

    private:
        /** @brief Add the edges of polygon a and b as hyperplanes of the polyhedron. The approximate normal is used when the polygon is degenerated */
        void addSourcePolygonEdgesAsHyperplanes(PluckerPolyhedron<P>* aPolyhedron, const GeometryConvexPolygon& aPolygon, bool clockWise, const MathVector3d& approximateNormal, VisibilityFeature::Kind aKind);

        /** @brief Add the extremal stabbing lines of a and b as hyperplanes of the polyhedron */
        void addExtremalStabbingLines(PluckerPolytope<P>* polytope, GeometryConvexPolygon& a, GeometryConvexPolygon& b, PluckerPolyhedron<P>* polyhedron);

        /** @brief Create the initial cell of a line family from a direction, the shadow of the first source clipped by the half planes of the second.

        The rays of the family that go through both sources are the rays whose screen
        point lies in the shadow of each of them, and their cell is the intersection of
        two convex polygons, a cell of two dimensions rather than four. No vertex pair
        of the sources is ever taken for a line here, which the general builder has to
        do: the direction of the light, and the vertices of the first source, are all
        the initialization asks for. */
        void addDirectionalScreenCell(PluckerPolytope<P>* polytope, GeometryConvexPolygon& a, GeometryConvexPolygon& b, PluckerPolyhedron<P>* polyhedron);

        /** @brief Create the edges from the facets representation to obtain the polytope skeleton*/
        void addEdges(PluckerPolytope<P>* polytope, PluckerPolyhedron<P>* aPolyhedron);

        MathPluckerNormalization<P, S> mNormalization;  /** < @brief How the points and the hyperplanes are normalized*/
        S mTolerance;         /** < @brief Tolerance used to compute the intersection with the Plucker quadric */
        VisibilityEventRegistry* mRegistry; /** < @brief Where the edges of the polygons are recorded, if any*/
        MathIntegerGrid mGrid;              /** < @brief The grid the polygons are snapped to: the lines through its points get exact predicates*/

        /** @brief Record in the polyhedron that the point anIndex is the line from a to b, when both are grid points*/
        void setInputLine(PluckerPolyhedron<P>* aPolyhedron, size_t anIndex, const MathVector3d& a, const MathVector3d& b)
        {
            if (mGrid.isOnGrid(a) && mGrid.isOnGrid(b))
            {
                aPolyhedron->setInputLine(anIndex, mGrid.toInteger(a), mGrid.toInteger(b));
            }
        }
    };

    template<class P, class S>
    inline PluckerPolytopeBuilder<P, S>::PluckerPolytopeBuilder(const MathPluckerNormalization<P, S>& aNormalization, S aTolerance, VisibilityEventRegistry* aRegistry, const MathIntegerGrid& aGrid)
        : mGrid(aGrid)
    {
        mNormalization = aNormalization;
        mTolerance = aTolerance;
        mRegistry = aRegistry;
    }

    template<class P, class S>
    inline void PluckerPolytopeBuilder<P, S>::addSourcePolygonEdgesAsHyperplanes(PluckerPolyhedron<P>* aPolyhedron, const GeometryConvexPolygon& aPolygon, bool clockWise, const MathVector3d& approximateNormal, VisibilityFeature::Kind aKind)
    {
        size_t n = aPolygon.getVertexCount();
        std::vector<MathVector3d> v;
        // which vertices are geometry: the padding of a degenerated polygon is not (v[1] always is, v[2] for a segment)
        std::vector<bool> isReal;

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
            isReal = { false, true, n == 2 };
        }
        else
        {
            // manage the case of a regular convex polygon
            v = aPolygon.getVertices();
            isReal.assign(v.size(), true);
        }
        size_t count = v.size();

        // Compute the oriented Plucker line representation of each edge and store it in Polyhedron
        for (size_t i = 0; i < count; i++)
        {
            const MathVector3d& v1 = v[i % count];
            const MathVector3d& v2 = v[(i + 1) % count];

            P  line = clockWise ? MathGeometry::createEdgeFromPoints<P>(v1, v2) : MathGeometry::createEdgeFromPoints<P>(v2, v1);
            line = mNormalization.normalizeHyperplane(line);
            const size_t index = aPolyhedron->add(line, ON_BOUNDARY, mNormalization.isUnitNorm(), mTolerance);
            if (clockWise)
            {
                setInputLine(aPolyhedron, index, v1, v2);
            }
            else
            {
                setInputLine(aPolyhedron, index, v2, v1);
            }

            if (mRegistry != nullptr)
            {
                VisibilityFeature feature;
                feature.mKind = aKind;
                feature.mEndpoints[0] = v1;
                feature.mEndpoints[1] = v2;
                feature.mRealEndpoint[0] = isReal[i % count];
                feature.mRealEndpoint[1] = isReal[(i + 1) % count];
                // the only real edge of a padded polygon joins its two real vertices
                feature.mSynthetic = !(feature.mRealEndpoint[0] && feature.mRealEndpoint[1]);
                mRegistry->registerHyperplane(index, feature);
            }
        }
    }

    template<class P, class S>
    inline void PluckerPolytopeBuilder<P, S>::addExtremalStabbingLines(PluckerPolytope<P> * polytope, GeometryConvexPolygon & a, GeometryConvexPolygon & b, PluckerPolyhedron<P> * aPolyhedron)
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
                myStabbingLine = mNormalization.normalizePoint(myStabbingLine);

                polytope->addExtremalStabbingLine(myStabbingLine);
                size_t vertexIndex = aPolyhedron->add(myStabbingLine, ON_BOUNDARY, mNormalization.isUnitNorm(), mTolerance);
                setInputLine(aPolyhedron, vertexIndex, v1, v2);

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

                aPolyhedron->initFacetsDescription(vertexIndex, myTemp);

                myStabbingNumber++;
            }
        }
    }

    template<class P, class S>
    inline PluckerPolytope<P>* PluckerPolytopeBuilder<P, S>::build(GeometryConvexPolygon & a, GeometryConvexPolygon & b, PluckerPolyhedron<P> * aPolyhedron)
    {
        PluckerPolytope<P>* myPolytope = new PluckerPolytope<P>();

        // The gravity centers are used to compute an approximated normal vector, that will be used when polygons are degenerated
        MathVector3d ga = MathGeometry::getGravityCenter(a);
        MathVector3d gb = MathGeometry::getGravityCenter(b);

        addSourcePolygonEdgesAsHyperplanes(aPolyhedron, a, true, gb - ga, VisibilityFeature::LIGHT_EDGE);
        addSourcePolygonEdgesAsHyperplanes(aPolyhedron, b, false, ga - gb, VisibilityFeature::RECEIVER_EDGE);

        if constexpr (MathPluckerCellTraits<P>::cAllLinesReal)
        {
            addDirectionalScreenCell(myPolytope, a, b, aPolyhedron);

            if (myPolytope->getVertices().empty())
            {
                /* The shadows do not meet: no ray of the light goes through both
                sources, and the empty cell is the answer, not an error waiting to
                be asserted. */
                return myPolytope;
            }
        }
        else
        {
            addExtremalStabbingLines(myPolytope, a, b, aPolyhedron);
            addEdges(myPolytope, aPolyhedron);
        }

        myPolytope->computeEdgesIntersectingQuadric(aPolyhedron, mTolerance);
        V_ASSERT(myPolytope->isValid(aPolyhedron, mNormalization.isUnitNorm(), mTolerance));

        return myPolytope;
    }

    template<class P, class S>
    inline void PluckerPolytopeBuilder<P, S>::addDirectionalScreenCell(PluckerPolytope<P>* polytope, GeometryConvexPolygon & a, GeometryConvexPolygon & b, PluckerPolyhedron<P> * aPolyhedron)
    {
        // The light runs from the center of the first source to the center of the second.
        MathVector3d ga = MathGeometry::getGravityCenter(a);
        MathVector3d gb = MathGeometry::getGravityCenter(b);
        MathVector3d d = gb - ga;
        V_ASSERT(!d.isZero(MathArithmetic<double>::Tolerance()));

        size_t countA = a.getVertexCount();
        size_t countB = b.getVertexCount();

        // The shadow of a degenerated source is a point or a line on the screen, and
        // the cell of such a light is the work of the days after this one.
        V_ASSERT(countA >= 3 && countB >= 3);

        const MathVector3_<S> myD = convert<MathVector3_<S>>(d);

        // The center ray of the second shadow lies strictly inside all its half planes,
        // and tells each of them which side the shadow is on.
        const P refB(myD, convert<MathVector3_<S>>(gb));

        /* The screen polygon starts as the shadow of the first source, the ray of the
        family through each of its vertices, and is clipped by the half planes of the
        edges of the second. The hyperplanes of both sources were added first, in order,
        by addSourcePolygonEdgesAsHyperplanes: the edges of b are planes countA to
        countA + countB. */
        std::vector<P> screen;
        for (size_t i = 0; i < countA; i++)
        {
            screen.push_back(P(myD, convert<MathVector3_<S>>(a.getVertex(i))));
        }

        for (size_t e = 0; e < countB && !screen.empty(); e++)
        {
            const P& plane = aPolyhedron->get(countA + e);

            const S ref = plane.dot(refB);
            V_ASSERT(MathArithmetic<S>::getAbs(ref) > mTolerance);
            const bool positive = ref > S(0.0);

            std::vector<P> clipped;
            const size_t n = screen.size();

            for (size_t i = 0; i < n; i++)
            {
                const P& cur = screen[i];
                const P& next = screen[(i + 1) % n];

                const S dc = plane.dot(cur);
                const S dn = plane.dot(next);

                const bool inCur = positive ? dc >= -mTolerance : dc <= mTolerance;
                const bool inNext = positive ? dn >= -mTolerance : dn <= mTolerance;

                if (inCur)
                {
                    clipped.push_back(cur);

                    if (!inNext)
                    {
                        clipped.push_back(MathGeometry::interpolate(dc, dn, cur, next, mTolerance));
                    }
                }
                else if (inNext)
                {
                    clipped.push_back(MathGeometry::interpolate(dc, dn, cur, next, mTolerance));
                }
            }

            screen.swap(clipped);
        }

        // The debris of a clip grazing an edge, two consecutive names of the same ray,
        // come down to one vertex.
        std::vector<P> cell;
        for (size_t i = 0; i < screen.size(); i++)
        {
            if (!cell.empty() && MathPredicates::isEdgeCollapsed(cell.back(), screen[i], mTolerance))
            {
                continue;
            }
            cell.push_back(screen[i]);
        }
        if (cell.size() > 1 && MathPredicates::isEdgeCollapsed(cell.front(), cell.back(), mTolerance))
        {
            cell.pop_back();
        }

        if (cell.size() < 3)
        {
            // Less than a polygon on the screen: the light has no cell to show.
            return;
        }

        /* A vertex of the cell lies on the boundary lines of the shadow it stands at,
        whichever source they come from: the incidence asks the Klein form of every
        hyperplane, and no bookkeeping of the clip is trusted with the facets. */
        const size_t planeCount = countA + countB;

        std::vector<size_t> indices;
        for (size_t i = 0; i < cell.size(); i++)
        {
            std::vector<size_t> facets;
            for (size_t f = 0; f < planeCount; f++)
            {
                if (MathPredicates::getRelativePosition(aPolyhedron->get(f).dot(cell[i]), mTolerance) == ON_BOUNDARY)
                {
                    facets.push_back(f);
                }
            }
            V_ASSERT(!facets.empty());

            size_t vertexIndex = aPolyhedron->add(cell[i], ON_BOUNDARY, mNormalization.isUnitNorm(), mTolerance);
            aPolyhedron->initFacetsDescription(vertexIndex, facets);

            // The vertices of the family are real rays themselves, the extremal
            // stabbing lines of the cell are its own corners.
            polytope->addExtremalStabbingLine(cell[i]);

            indices.push_back(vertexIndex);
        }

        // The edges join consecutive corners of the screen polygon, each pair of them
        // held by the boundary line they lie on together.
        for (size_t i = 0; i < indices.size(); i++)
        {
            polytope->addEdge(indices[i], indices[(i + 1) % indices.size()], aPolyhedron);
        }
    }

    template<class P, class S>
    inline void PluckerPolytopeBuilder<P, S>::addEdges(PluckerPolytope<P> * polytope, PluckerPolyhedron<P> * aPolyhedron)
    {
        // We iterate through all the extremal stabbing lines, and creates and edge for the vertices that share at least three common facets
        // The complexity is O(n^2), and more optimized version could be devised but we keep this one for simplicity
        for (int i = 0; i < aPolyhedron->getLinesCount(); i++)
        {
            for (int j = i + 1; j < aPolyhedron->getLinesCount(); j++)
            {
                if (MathCombinatorial::haveAtLeastNCommonFacets(aPolyhedron->getFacetsDescription(i), aPolyhedron->getFacetsDescription(j), MathPluckerCellTraits<P>::cCommonFacetsForEdge))
                {
                    polytope->addEdge(i, j, aPolyhedron);
                }
            }
        }
    }
}
