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

#include <string>
#include <iostream>

/* The umbrella first, as every consumer of the library takes it: the headers of the
space are caught in a circle of includes otherwise. */
#include "visilib.h"

#include "math_combinatorial.h"

#include "math_plucker_parallel_6.h"
#include "plucker_polytope_builder.h"
#include "plucker_polytope_splitter.h"

#ifdef ENABLE_GMP
#include "math_gmp_rational.h"
#include "math_gmp_float.h"
#endif
#ifdef ENABLE_MPFR
#include "math_mpfr.h"
#endif

using namespace visilib;

bool MathCombinatorialTest(std::string& )
{
	{//MathCombinatorial::hasFacet
		if (MathCombinatorial::hasFacet(std::vector<size_t>(), 1))
			{ std::cout << "Error in line " <<  __LINE__ -1 << std::endl; return false;}

		if (MathCombinatorial::hasFacet(std::vector<size_t>{ 1, 2, 3, 4 }, 0))
			{ std::cout << "Error in line " <<  __LINE__ -1 << std::endl; return false;}

		for (size_t i = 1; i <= 6; i++)
		{
			if (!MathCombinatorial::hasFacet(std::vector<size_t>{ 1, 2, 3, 4, 5, 6 }, i))
				{ std::cout << "Error in line " <<  __LINE__ -1 << std::endl; return false;}
		}
	}

	{//MathCombinatorial::haveAtLeastNCommonFacets
		for (size_t i = 0; i < 5; i++)
		{
			if (MathCombinatorial::haveAtLeastNCommonFacets(std::vector<size_t>{1, 2, 3}, std::vector<size_t>{1, 2, 3}, i) != (i <= 3))
				{ std::cout << "Error in line " <<  __LINE__ -1 << std::endl; return false;}

			if (MathCombinatorial::haveAtLeastNCommonFacets(std::vector<size_t>{1, 2, 3, 4, 5, 6}, std::vector<size_t>{1, 2, 3}, i) != (i <= 3))
				{ std::cout << "Error in line " <<  __LINE__ -1 << std::endl; return false;}

			if (MathCombinatorial::haveAtLeastNCommonFacets(std::vector<size_t>{1, 2, 3}, std::vector<size_t>{4, 5, 6, 7, 8, 9 }, i) != (i <= 0))
				{ std::cout << "Error in line " <<  __LINE__ -1 << std::endl; return false;}

			if (MathCombinatorial::haveAtLeastNCommonFacets(std::vector<size_t>{1, 6, 9, 12, 15, 24, 31}, std::vector<size_t>{6, 8, 10, 22, 24, 25, 31, 32, 33 }, i) != (i <= 3))
				{ std::cout << "Error in line " <<  __LINE__ -1 << std::endl; return false;}

			if (MathCombinatorial::haveAtLeastNCommonFacets(std::vector<size_t>{1, 12, 13, 14, 15, 16}, std::vector<size_t>{1, 2, 3, 4, 5, 6, 7, 16  }, i) != (i <= 2))
				{ std::cout << "Error in line " <<  __LINE__ -1 << std::endl; return false;}

			if (MathCombinatorial::haveAtLeastNCommonFacets(std::vector<size_t>{5, 10, 15, 20, 25, 30}, std::vector<size_t>{11, 12, 13, 14, 15, 16, 17, 18,19  }, i) != (i <= 1))
				{ std::cout << "Error in line " <<  __LINE__ -1 << std::endl; return false;}

			if (MathCombinatorial::haveAtLeastNCommonFacets(std::vector<size_t>{1, 1, 1, 1}, std::vector<size_t>{1, 1, 1, 1, 1, 1 }, i) != (i <= 4))
				{ std::cout << "Error in line " <<  __LINE__ -1 << std::endl; return false;}

			if (MathCombinatorial::haveAtLeastNCommonFacets(std::vector<size_t>{}, std::vector<size_t>{1, 2, 3, 4 }, i) != (i <= 0))
				{ std::cout << "Error in line " <<  __LINE__ -1 << std::endl; return false;}

			if (MathCombinatorial::haveAtLeastNCommonFacets(std::vector<size_t>{1, 2, 3, 4}, std::vector<size_t>{}, i) != (i <= 0))
				{ std::cout << "Error in line " <<  __LINE__ -1 << std::endl; return false;}

			if (MathCombinatorial::haveAtLeastNCommonFacets(std::vector<size_t>{}, std::vector<size_t>{}, i) != (i <= 0))
				{ std::cout << "Error in line " <<  __LINE__ -1 << std::endl; return false;}
		}
	}

	{//MathCombinatorial::initFacets
		std::vector<size_t> result;
		MathCombinatorial::initFacets(std::vector<size_t>{1, 2, 3}, std::vector<size_t>{1, 2, 3}, 4, result);
		if (result != std::vector<size_t>{ 1, 2, 3, 4 })
			{ std::cout << "Error in line " <<  __LINE__ -1 << std::endl; return false;}

		MathCombinatorial::initFacets(std::vector<size_t>{1, 2, 3}, std::vector<size_t>{1, 2, 3, 4, 5, 6, 7}, 8, result);
		if (result != std::vector<size_t>{ 1, 2, 3, 8 })
			{ std::cout << "Error in line " <<  __LINE__ -1 << std::endl; return false;}

		MathCombinatorial::initFacets(std::vector<size_t>{1, 2, 3, 4, 5, 6, 7}, std::vector<size_t>{1, 2, 3}, 8, result);
		if (result != std::vector<size_t>{ 1, 2, 3, 8 })
			{ std::cout << "Error in line " <<  __LINE__ -1 << std::endl; return false;}

		MathCombinatorial::initFacets(std::vector<size_t>{1, 2, 3, 4, 5, 6}, std::vector<size_t>{2, 4, 6, 8, 10}, 11, result);
		if (result != std::vector<size_t>{ 2, 4, 6, 11 })
			{ std::cout << "Error in line " <<  __LINE__ -1 << std::endl; return false;}

		MathCombinatorial::initFacets(std::vector<size_t>{1, 2, 3, 4, 5, 6}, std::vector<size_t>{1, 2, 3, 4, 5, 6}, 7, result);
		if (result != std::vector<size_t>{ 1, 2, 3, 4, 5, 6, 7 })
			{ std::cout << "Error in line " <<  __LINE__ -1 << std::endl; return false;}
	}


    std::cout << "MathCombinatorialTest SUCCESS" << std::endl;

	return true;
}

namespace
{
    /* The vectors of the checks are all of integer, or half, coordinates: the moments,
    the projections and the Klein products stay exact whatever the arithmetic, so the
    same expected values hold for double, float and the exact types alike. The square
    rooted checks, normalization and distance, ask of the wrapper arithmetics what
    they do not promise, and are left to the floating ones. */

    template<class S>
    MathVector3_<S> v3(double x, double y, double z)
    {
        return convert<MathVector3_<S>>(MathVector3d(x, y, z));
    }

    bool fail(const int line)
    {
        std::cout << "Error in line " << line << std::endl;
        return false;
    }

    template<class S>
    bool isNought(const S& value, const S& tolerance)
    {
        return MathArithmetic<S>::getAbs(value) <= tolerance;
    }

    /* The exact arithmetics compare with < and <= only, so equality is asked as the
    vanishing of a difference; the data of the checks being exact, a nought
    difference is a nought difference whatever the model. */
    template<class S>
    bool sameVector(const MathVector3_<S>& a, const MathVector3_<S>& b, const S& tolerance)
    {
        return MathVector3_<S>(a - b).isZero(tolerance);
    }

    template<class P, class S>
    bool checkParallelAlgebra(bool checkSqrt)
    {
        const S tolerance = MathArithmetic<S>::Tolerance();

        const MathVector3_<S> dx = v3<S>(1, 0, 0);

        // The rays of the family, along x, through three points.
        const P ray1(dx, v3<S>(2, 3, 5));
        const P ray2(dx, v3<S>(2, -1, 4));
        const P ray3(dx, v3<S>(2, 7, 0));

        // The moment is the cross product, the screen point recovers the projected point,
        // and the quadric is identically nought.
        if (!sameVector(ray1.getLocation(), v3<S>(0, 5, -3), tolerance)) return fail(__LINE__);
        if (!sameVector(ray1.getScreenPoint(), v3<S>(0, 3, 5), tolerance)) return fail(__LINE__);
        if (!isNought(ray1.getQuadric(), tolerance)) return fail(__LINE__);

        // An edge carried by the class: its plane against the rays of the family.
        const P edge = P::throughEdge(v3<S>(0, 1, 1), v3<S>(0, -1, 1));
        const P skewRay(dx, v3<S>(0, 0, 2));
        const P meetRay(dx, v3<S>(0, 0, 1));

        if (!isNought(edge.getQuadric(), tolerance)) return fail(__LINE__);
        if (!isNought(edge.dot(edge), tolerance)) return fail(__LINE__);
        if (!isNought(edge.dot(meetRay), tolerance)) return fail(__LINE__);
        if (!isNought(edge.dot(skewRay) + S(2.0), tolerance)) return fail(__LINE__);

        // The accumulation of the rays from the zero, the gravity center of the family
        // on its screen: what the polytopes are made of, with no correction to chase.
        P acc = P::Zero();
        if (!acc.isZero(tolerance)) return fail(__LINE__);
        acc += ray1;
        acc += ray2;
        acc += ray3;
        acc = acc / S(3.0);

        if (!sameVector(acc.getDirection(), dx, tolerance)) return fail(__LINE__);
        if (!sameVector(acc.getLocation(), v3<S>(0, 3, -3), tolerance)) return fail(__LINE__);
        if (!sameVector(acc.getScreenPoint(), v3<S>(0, 3, 3), tolerance)) return fail(__LINE__);
        if (!isNought(acc.getQuadric(), tolerance)) return fail(__LINE__);

        // A sum across families, where the projection of the moment is what keeps the
        // result a line at all.
        const P rayY(v3<S>(0, 1, 0), v3<S>(1, 0, 0));
        const P mixed = ray1 + rayY;
        if (!sameVector(mixed.getDirection(), v3<S>(1, 1, 0), tolerance)) return fail(__LINE__);
        if (!sameVector(mixed.getLocation(), v3<S>(-2.5, 2.5, -2), tolerance)) return fail(__LINE__);
        if (!isNought(mixed.getQuadric(), tolerance)) return fail(__LINE__);

        // Scaling both parts renormalizes the coordinates and does not move the ray.
        const P doubled = ray1 * S(2.0);
        if (!sameVector(doubled.getLocation(), v3<S>(0, 10, -6), tolerance)) return fail(__LINE__);
        if (!sameVector(doubled.getScreenPoint(), ray1.getScreenPoint(), tolerance)) return fail(__LINE__);

        // The gradient of the quadric swaps the parts, and satisfies the same invariant.
        const P grad = ray1.quadricGradient();
        if (!sameVector(grad.getDirection(), ray1.getLocation(), tolerance)) return fail(__LINE__);
        if (!sameVector(grad.getLocation(), ray1.getDirection(), tolerance)) return fail(__LINE__);
        if (!isNought(grad.getQuadric(), tolerance)) return fail(__LINE__);

        // The ordinary six dimensional norm: 1 + 25 + 9.
        if (!isNought(ray1.scalarDot(ray1) - S(35.0), tolerance)) return fail(__LINE__);

        // A copy, and the output operator the assertions read through.
        const P copy(ray1);
        if (!sameVector(copy.getLocation(), ray1.getLocation(), tolerance)) return fail(__LINE__);
        if (!sameVector(copy.getDirection(), ray1.getDirection(), tolerance)) return fail(__LINE__);

        if (checkSqrt)
        {
            const P one = ray1.getNormalized();
            if (!isNought(one.scalarDot(one) - S(1.0), tolerance * 100)) return fail(__LINE__);
            if (!(one.getScreenPoint() - ray1.getScreenPoint()).isZero(tolerance * 100)) return fail(__LINE__);

            // The offset of two rays of the family, their screen points one unit apart.
            const P next(dx, v3<S>(0, 3, 6));
            if (!isNought(ray1.getDistance(next) - S(1.0), tolerance * 100)) return fail(__LINE__);
        }

        return true;
    }

    /* The supporting types of the exact query, taken with the family as their point
    type. The pipeline itself is not run: what is checked here is that the machinery
    of the space accepts the class, and keeps its predicates and its containers honest
    on the values of the family. */
    bool checkParallelQueryTypes()
    {
        typedef MathPluckerParallel6d ParallelLine;
        typedef PluckerPolyhedron<ParallelLine> ParallelPolyhedron;
        typedef PluckerPolytope<ParallelLine> ParallelPolytope;
        typedef VisibilityExactQuery_<ParallelLine, double> ParallelQuery;

        const double tolerance = MathArithmetic<double>::Tolerance();

        ParallelPolyhedron polyhedron;
        const ParallelLine edge = ParallelLine::throughEdge(MathVector3d(0, 1, 1), MathVector3d(0, -1, 1));
        const ParallelLine ray(MathVector3d(1, 0, 0), MathVector3d(0, 0, 2));

        const size_t index = polyhedron.add(edge, ON_BOUNDARY, false, tolerance);
        if (index != 0) return fail(__LINE__);
        if (polyhedron.getLinesCount() != 1) return fail(__LINE__);
        if (polyhedron.getQuadricRelativePosition(0) != ON_BOUNDARY) return fail(__LINE__);
        if (polyhedron.get(0).dot(ray) >= -tolerance) return fail(__LINE__);

        ParallelPolytope polytope;
        polytope.setRepresentativeLine(ray);
        polytope.addExtremalStabbingLine(ray);
        if (polytope.getExtremalStabbingLinesCount() != 1) return fail(__LINE__);
        if (polytope.getExtremalStabbingLine(0).getLocation() != ray.getLocation()) return fail(__LINE__);

        // The predicates of the space, read on the planes and the rays of the family.
        if (MathPredicates::getVertexPlaneRelativePosition(edge, ray, tolerance) != ON_NEGATIVE_SIDE) return fail(__LINE__);
        if (MathPredicates::getQuadricRelativePosition(ray, tolerance) != ON_BOUNDARY) return fail(__LINE__);

        // The query takes the family as its point type: the definition stands, the run
        // of the pipeline is the work of the days after this one.
        (void)sizeof(ParallelQuery);

        return true;
    }

    /* The directional initialization: the cell of the light through two sources,
    built from the direction the centers tell and the vertices of the first source,
    clipped by the half planes of the second. No pair of vertices is ever taken for
    a line here, which is what the general builder does. */
    bool checkDirectionalScreenCell()
    {
        typedef MathPluckerParallel6d ParallelLine;
        typedef PluckerPolyhedron<ParallelLine> ParallelPolyhedron;
        typedef PluckerPolytope<ParallelLine> ParallelPolytope;

        const double tolerance = MathArithmetic<double>::Tolerance();

        // Two unit squares facing each other along x, centers (0,0,0) and (2,0,0).
        const float squareA[] = { 0, -1, -1, 0, 1, -1, 0, 1, 1, 0, -1, 1 };
        const float squareB[] = { 2, -1, -1, 2, 1, -1, 2, 1, 1, 2, -1, 1 };

        // The same square, half the size: its shadow lies within the first one.
        const float squareS[] = { 2, -0.5, -0.5, 2, 0.5, -0.5, 2, 0.5, 0.5, 2, -0.5, 0.5 };
        const float squareW[] = { 0, -2, -2, 0, 2, -2, 0, 2, 2, 0, -2, 2 };

        // Two bars crossing at right angles, their shadows overlapping in a square.
        const float barA[] = { 0, -2, -0.5, 0, 2, -0.5, 0, 2, 0.5, 0, -2, 0.5 };
        const float barB[] = { 2, -0.5, -2, 2, 0.5, -2, 2, 0.5, 2, 2, -0.5, 2 };

        // Coincident shadows: the cell is the square, four corners each lying on two
        // lines of either source, and the center ray answers for it.
        {
            GeometryConvexPolygon a(squareA, 4);
            GeometryConvexPolygon b(squareB, 4);

            ParallelPolyhedron polyhedron;
            PluckerPolytopeBuilder<ParallelLine, double> builder(false, tolerance);
            ParallelPolytope* polytope = builder.build(a, b, &polyhedron);

            if (polytope->getVertices().size() != 4) return fail(__LINE__);
            if (polytope->getEdgeCount() != 4) return fail(__LINE__);
            if (polytope->getExtremalStabbingLinesCount() != 4) return fail(__LINE__);
            if (!polytope->containsRealLines()) return fail(__LINE__);

            const RepresentativeLine<ParallelLine> rep = MathGeometry::representativeLine<ParallelLine, double>(polytope, &polyhedron, tolerance, true);
            if (!rep.insideCell) return fail(__LINE__);
            if (!sameVector(rep.line.getScreenPoint(), MathVector3d(0, 0, 0), tolerance)) return fail(__LINE__);

            delete polytope;
        }

        // A shadow within the other: the cell is the small square, four corners.
        {
            GeometryConvexPolygon a(squareW, 4);
            GeometryConvexPolygon b(squareS, 4);

            ParallelPolyhedron polyhedron;
            PluckerPolytopeBuilder<ParallelLine, double> builder(false, tolerance);
            ParallelPolytope* polytope = builder.build(a, b, &polyhedron);

            if (polytope->getVertices().size() != 4) return fail(__LINE__);
            if (polytope->getEdgeCount() != 4) return fail(__LINE__);

            delete polytope;
        }

        /* The crossing bars: the cell is the square of the overlap, its four corners
        each held by a line of either source. A silhouette plane then cuts it, the
        very split of the solver, and both children keep a certified representative. */
        {
            GeometryConvexPolygon a(barA, 4);
            GeometryConvexPolygon b(barB, 4);

            ParallelPolyhedron polyhedron;
            PluckerPolytopeBuilder<ParallelLine, double> builder(false, tolerance);
            ParallelPolytope* polytope = builder.build(a, b, &polyhedron);

            if (polytope->getVertices().size() != 4) return fail(__LINE__);
            if (polytope->getEdgeCount() != 4) return fail(__LINE__);

            const RepresentativeLine<ParallelLine> rep = MathGeometry::representativeLine<ParallelLine, double>(polytope, &polyhedron, tolerance, true);
            if (!rep.insideCell) return fail(__LINE__);
            if (!sameVector(rep.line.getScreenPoint(), MathVector3d(0, 0, 0), tolerance)) return fail(__LINE__);

            // A screen line at y = 0.25, carried by an edge standing in the light.
            const ParallelLine silPlane = ParallelLine::throughEdge(MathVector3d(0, 0.25, -5), MathVector3d(0, 0.25, 5));
            const size_t silIndex = polyhedron.add(silPlane, ON_BOUNDARY, false, tolerance);

            ParallelPolytope* left = new ParallelPolytope();
            ParallelPolytope* right = new ParallelPolytope();

            const GeometryPositionType position = PluckerPolytopeSplitter<ParallelLine, double>::split(&polyhedron, silPlane, polytope, left, right, silIndex, false, tolerance);
            if (position != ON_BOTH_SIDES) return fail(__LINE__);

            left->computeEdgesIntersectingQuadric(&polyhedron, tolerance);
            right->computeEdgesIntersectingQuadric(&polyhedron, tolerance);

            if (!left->containsRealLines()) return fail(__LINE__);
            if (!right->containsRealLines()) return fail(__LINE__);
            if (left->getVertices().size() != 4) return fail(__LINE__);
            if (right->getVertices().size() != 4) return fail(__LINE__);

            const RepresentativeLine<ParallelLine> repLeft = MathGeometry::representativeLine<ParallelLine, double>(left, &polyhedron, tolerance, true);
            const RepresentativeLine<ParallelLine> repRight = MathGeometry::representativeLine<ParallelLine, double>(right, &polyhedron, tolerance, true);
            if (!repLeft.insideCell) return fail(__LINE__);
            if (!repRight.insideCell) return fail(__LINE__);

            // The children stand on either side of the plane that split them.
            if (MathPredicates::getRelativePosition(left, &polyhedron, silPlane, tolerance) != ON_NEGATIVE_SIDE) return fail(__LINE__);
            if (MathPredicates::getRelativePosition(right, &polyhedron, silPlane, tolerance) != ON_POSITIVE_SIDE) return fail(__LINE__);

            delete polytope;
            delete left;
            delete right;
        }

        return true;
    }
}

bool MathPluckerParallelTest(std::string& )
{
    if (!checkParallelAlgebra<MathPluckerParallel6f, float>(true))
        return false;
    if (!checkParallelAlgebra<MathPluckerParallel6d, double>(true))
        return false;
#ifdef ENABLE_GMP
    if (!checkParallelAlgebra<MathPluckerParallel6<MathGmpRational>, MathGmpRational>(false))
        return false;
    if (!checkParallelAlgebra<MathPluckerParallel6<MathGmpFloat>, MathGmpFloat>(false))
        return false;
#endif
#ifdef ENABLE_LEDA
    if (!checkParallelAlgebra<MathPluckerParallel6<MathLedaReal>, MathLedaReal>(false))
        return false;
#endif
#ifdef ENABLE_MPFR
    MathMpfr::set_default_precision(MathArithmetic<MathMpfr>::bitsCount());
    if (!checkParallelAlgebra<MathPluckerParallel6<MathMpfr>, MathMpfr>(false))
        return false;
#endif

    if (!checkParallelQueryTypes())
        return false;
    if (!checkDirectionalScreenCell())
        return false;

    std::cout << "MathPluckerParallelTest SUCCESS" << std::endl;

    return true;
}
