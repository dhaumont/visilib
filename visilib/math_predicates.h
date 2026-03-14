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

#include "math_plucker_6.h"
#include "geometry_position_type.h"

namespace visilib
{
    class GeometryConvexPolygon;
    template<class S>
    class PluckerPolytope;

    template<class S>
    class PluckerPolyhedron;

    /** @brief Contains geometry predicates function */

    class MathPredicates
    {
    public:
        template <class P>          static bool                 isPointOnQuadric(const P& point);
        template <class P, class S> static bool                 isEdgeCollapsed(const P& a, const P& b, S tolerance);
        template <class P, class S> static bool                 isNormalized(const P& point, S tolerance);
        template <class S>          static GeometryPositionType getRelativePosition(S dot, S epsilon);
        template <class P, class S> static GeometryPositionType getVertexPlaneRelativePosition(const P& plane, const P& point, S epsilon);
        template <class P, class S> static GeometryPositionType getQuadricRelativePosition(const P& point, S epsilon);
        template <class P>          static bool                 hasPluckerPolytopeIntersectionWithQuadric(PluckerPolytope<P>* polytope, PluckerPolyhedron<P>* polyhedron);
        template <class P>          static GeometryPositionType getRelativePosition(PluckerPolytope<P>* polytope, PluckerPolyhedron<P>* polyhedron, const P& aPlane0, const P& aPlane1, const P& aPlane2);
        template <class P, class S> static GeometryPositionType getRelativePosition(PluckerPolytope<P>* polytope, PluckerPolyhedron<P>* polyhedron, const P& aPlane, S tolerance);

        static GeometryPositionType getRelativePosition(const std::vector<MathVector3d>& points, const MathPlane3d& aPlane);
    };

    template <class P>
    inline bool MathPredicates::isPointOnQuadric(const P& a)
    {
        return getRelativePosition(a, a) == ON_BOUNDARY;
    }

    template <class P, class S>
    inline bool MathPredicates::isEdgeCollapsed(const P& a, const P& b, S tolerance)
    {
      P difference(a);
      difference -= b;
      return difference.isZero(tolerance);
    }

    template <class P, class S>
    inline GeometryPositionType MathPredicates::getVertexPlaneRelativePosition(const P& plane, const P& point, S tolerance)
    {
        auto result = plane.dot(point);
        return getRelativePosition(result, tolerance);
    }

    template <class S>
    inline GeometryPositionType MathPredicates::getRelativePosition(S dot, S epsilon)
    {        
        if (dot < -epsilon)
        {
            return ON_NEGATIVE_SIDE;
        }
        else if (dot > epsilon)
        {
            return ON_POSITIVE_SIDE;
        }

        return ON_BOUNDARY;
    }
    template <class P, class S>
    inline GeometryPositionType MathPredicates::getQuadricRelativePosition(const P& a, S tolerance)
    {
        return getVertexPlaneRelativePosition(a, a, tolerance);
    }

    template <class P>
    inline bool MathPredicates::hasPluckerPolytopeIntersectionWithQuadric(PluckerPolytope<P>* polytope, PluckerPolyhedron<P>* polyhedron)
    {
        bool hasPointOnTheLeft = false;
        bool hasPointOnTheRight = false;

        auto myVertices = polytope->getVertices();

        for (auto iter = myVertices.begin(); iter != myVertices.end(); iter++)
        {
            size_t v = *iter;
            GeometryPositionType position = polyhedron->getQuadricRelativePosition(v);

            if (position == ON_BOUNDARY)
            {
                return true;
            }
            else if (position == ON_NEGATIVE_SIDE)
            {
                hasPointOnTheLeft = true;
            }
            else
            {
                V_ASSERT(position == ON_POSITIVE_SIDE);
                hasPointOnTheRight = true;
            }

            if (hasPointOnTheLeft && hasPointOnTheRight)
            {
                return true;
            }
        }

        V_ASSERT(!hasPointOnTheLeft || !hasPointOnTheRight);

        if (hasPointOnTheLeft)
        {
            return false;
        }
        else if (hasPointOnTheRight)
        {
            return  false;
        }

        return true;
    }

    template<typename P, class S>
    inline bool MathPredicates::isNormalized(const P & a, S tolerance)
    {
        //P norm = a.R6DotProduct(a);
        auto norm = a.getLocation().getSquaredNorm()+a.getDirection().getSquaredNorm();
        return (MathArithmetic<S>::getAbs(norm - (S)1.0) <= tolerance);
    }

    inline GeometryPositionType MathPredicates::getRelativePosition(const std::vector<MathVector3d> & points, const MathPlane3d & aPlane)
    {
        int	myRight = 0;
        int	myLeft = 0;

        for (size_t i = 0; i < points.size(); i++)
        {
            double d = aPlane.dot(points[i]);

            if (d > MathArithmetic<double>::Tolerance())
            {
                myRight++;
            }
            else if (d < -MathArithmetic<double>::Tolerance())
            {
                myLeft++;
            }
            if (myRight > 0 && myLeft > 0)
            {
                return ON_BOUNDARY;
            }
        }

        if (myRight == 0 && myLeft == 0 && points.size() > 0)
        {
            // Relative position is on the boundary if all points were on the plane.
            return ON_BOUNDARY;
        }
        else if (myRight == 0)
        {
            V_ASSERT(myLeft > 0);
            return ON_NEGATIVE_SIDE;
        }
        else if (myLeft == 0)
        {
            V_ASSERT(myRight > 0);
            return ON_POSITIVE_SIDE;
        }
        V_ASSERT(0);
        return ON_BOUNDARY;
    }

    template <class P>
    inline GeometryPositionType MathPredicates::getRelativePosition(PluckerPolytope<P> * polytope, PluckerPolyhedron<P> * polyhedron, const P & aPlane0, const P & aPlane1, const P & aPlane2)
    {
        GeometryPositionType position0 = getRelativePosition(polytope, polyhedron, aPlane0);
        if (position0 == ON_POSITIVE_SIDE)
            return ON_POSITIVE_SIDE;
        GeometryPositionType position1 = getRelativePosition(polytope, polyhedron, aPlane1);
        if (position1 == ON_POSITIVE_SIDE)
            return ON_POSITIVE_SIDE;
        GeometryPositionType position2 = getRelativePosition(polytope, polyhedron, aPlane2);
        if (position2 == ON_POSITIVE_SIDE)
            return ON_POSITIVE_SIDE;

        if (position0 == ON_NEGATIVE_SIDE && position1 == ON_NEGATIVE_SIDE && position2 == ON_NEGATIVE_SIDE)
            return ON_NEGATIVE_SIDE;

        return ON_BOUNDARY;
    }

    template <class P, class S>
    inline GeometryPositionType MathPredicates::getRelativePosition(PluckerPolytope<P> * aPolytope, PluckerPolyhedron<P> * aPolyhedron, const P & aPlane, S tolerance)
    {
        bool hasPointOnTheLeft = false;
        bool hasPointOnTheRight = false;
        bool hasPointOnBoundary = false;

        auto myVertices = aPolytope->getVertices();

        for (auto iter = myVertices.begin(); iter != myVertices.end(); iter++)
        {
            int v = *iter;
            
            S result = aPlane.dot(aPolyhedron->get(v));
            GeometryPositionType position = MathPredicates::getRelativePosition(result, tolerance);
            V_ASSERT(position == MathPredicates::getVertexPlaneRelativePosition(aPlane, aPolyhedron->get(v), tolerance));
           
            if (position == ON_NEGATIVE_SIDE)
            {
                hasPointOnTheLeft = true;
            }
            else if (position == ON_POSITIVE_SIDE)
            {
                V_ASSERT(position == ON_POSITIVE_SIDE);
                hasPointOnTheRight = true;
            }
            else
            {
                V_ASSERT(position == ON_BOUNDARY);
                hasPointOnBoundary = true;            
            }

            if (hasPointOnTheLeft && hasPointOnTheRight)
            {
                return ON_BOTH_SIDES;
            }
        }

        V_ASSERT(!hasPointOnTheLeft || !hasPointOnTheRight);

        if (!hasPointOnTheLeft && !hasPointOnTheRight)
        {
            V_ASSERT(hasPointOnBoundary);
            return ON_BOUNDARY;
        }
        if (!hasPointOnTheLeft)
        {
            return ON_POSITIVE_SIDE;
        }
        else if (!hasPointOnTheRight)
        {
            return  ON_NEGATIVE_SIDE;
        }
        V_ASSERT(0);
        return ON_BOTH_SIDES;
    }
}
