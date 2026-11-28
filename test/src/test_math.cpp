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
#include <cmath>

/* The umbrella first, as every consumer of the library takes it: the headers of the
space are caught in a circle of includes otherwise. */
#include "visilib.h"

#include "math_combinatorial.h"
#include "math_conic_2.h"
#include "math_frame_3.h"
#include "math_polynomial.h"
#include "math_integer_predicates.h"
#include "math_plucker_normalization.h"
#include "math_sos.h"

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
        RepresentativeLine<ParallelLine> representative;
        representative.line = ray;
        polytope.setCachedRepresentativeLine(&polyhedron, true, representative);
        if (polytope.getRepresentativeLine().getLocation() != ray.getLocation()) return fail(__LINE__);
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
            PluckerPolytopeBuilder<ParallelLine, double> builder(MathPluckerNormalization<ParallelLine, double>(), tolerance);
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
            PluckerPolytopeBuilder<ParallelLine, double> builder(MathPluckerNormalization<ParallelLine, double>(), tolerance);
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
            PluckerPolytopeBuilder<ParallelLine, double> builder(MathPluckerNormalization<ParallelLine, double>(), tolerance);
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

            const GeometryPositionType position = PluckerPolytopeSplitter<ParallelLine, double>::split(&polyhedron, silPlane, polytope, left, right, silIndex, MathPluckerNormalization<ParallelLine, double>(), tolerance);
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

    /* The beam capsule: the convex hull of two balls, one at each end of the axis with the
    radius of that end, which the aperture finder casts wide when its thin line misses — the
    body is made to hold every line of a cell, so that casting it can miss no occluder.

    All the hit and miss points below are dyadic: the float and the double runs compute the
    very same numbers exactly, and no contact asked of the test sits closer to its grazing
    than a sixteenth of a radius, which neither arithmetic rounds across. The nearest-point
    checks use one tenth and six tenths, and are held to a tolerance rather than to nought. */
    template <class S>
    static bool checkBeamCapsule()
    {
        typedef MathVector3_<S> V;

        const V a = v3<S>(0, 0, 0);
        const V b = v3<S>(0, 0, 1);

        auto hits = [](const V& start, const V& end, S r0, S r1, const V& t0, const V& t1, const V& t2)
        {
            return MathGeometry::hitsBeamCapsule<S>(start, end, r0, r1, t0, t1, t2);
        };

        const S r0(0.125);
        const S r1(0.25);

        /* The crossing face: the beam measures 0.1875 at mid-height. */
        {
            const V t0 = v3<S>(0.125, 0, 0.5);
            const V t1 = v3<S>(0.5, 0.25, 0.5);
            const V t2 = v3<S>(0.5, -0.25, 0.5);
            if (!hits(a, b, r0, r1, t0, t1, t2)) return fail(__LINE__);
        }
        {
            const V t0 = v3<S>(0.25, 0, 0.5);
            const V t1 = v3<S>(0.5, 0.25, 0.5);
            const V t2 = v3<S>(0.5, -0.25, 0.5);
            if (hits(a, b, r0, r1, t0, t1, t2)) return fail(__LINE__);
        }

        /* The face lying in a plane parallel to the axis: the old test of the cylinder pierced
        the plane with the axis and called such a face a miss however close it stood — a
        grazing occluder, which is the one a shadow is made of. The ball answers it plainly:
        the axis is 0.125 from the face, the beam is wider there. */
        {
            const V t0 = v3<S>(0.125, 0, 0.4375);
            const V t1 = v3<S>(0.125, 0, 0.5625);
            const V t2 = v3<S>(0.5, 0, 0.5);
            if (!hits(a, b, r0, r1, t0, t1, t2)) return fail(__LINE__);
        }

        /* Beyond the caps: the beam is a body of measure, not an infinite cylinder — a face
        lying past the wide end, four tenths of an axis beyond it, stands out of every ball. */
        {
            const V t0 = v3<S>(0.125, 0, 1.4);
            const V t1 = v3<S>(-0.125, 0.125, 1.4);
            const V t2 = v3<S>(-0.125, -0.125, 1.4);
            if (hits(a, b, r0, r1, t0, t1, t2)) return fail(__LINE__);
        }

        /* Through the cap of the start, on the far side of the axis and yet inside the ball. */
        {
            const V t0 = v3<S>(0.0625, 0, -0.0625);
            const V t1 = v3<S>(-0.125, 0.125, -0.1875);
            const V t2 = v3<S>(-0.125, -0.125, -0.1875);
            if (!hits(a, b, r0, r1, t0, t1, t2)) return fail(__LINE__);
        }

        /* Near the wide end: at the cap itself the beam is a quarter wide. */
        {
            const V t0 = v3<S>(0.2, 0, 1);
            const V t1 = v3<S>(0.5, 0.25, 1);
            const V t2 = v3<S>(0.5, -0.25, 1);
            if (!hits(a, b, r0, r1, t0, t1, t2)) return fail(__LINE__);
        }
        {
            const V t0 = v3<S>(0.3, 0, 1);
            const V t1 = v3<S>(0.5, 0.25, 1);
            const V t2 = v3<S>(0.5, -0.25, 1);
            if (hits(a, b, r0, r1, t0, t1, t2)) return fail(__LINE__);
        }

        /* A beam drawn to nothing is the axis itself: a face straddling the segment is held,
        a face beside it is not. */
        {
            const V t0 = v3<S>(0.0625, 0.0625, 0.5);
            const V t1 = v3<S>(-0.0625, 0.0625, 0.5);
            const V t2 = v3<S>(0, -0.0625, 0.5);
            if (!hits(a, b, S(0), S(0), t0, t1, t2)) return fail(__LINE__);

            const V g0 = v3<S>(0.0625, 0.0625, 0.5);
            const V g1 = v3<S>(0.5, 0.0625, 0.5);
            const V g2 = v3<S>(0.5, -0.0625, 0.5);
            if (hits(a, b, S(0), S(0), g0, g1, g2)) return fail(__LINE__);
        }

        /* The degenerate axis: the beam is the bigger of the two balls it wanted to be. */
        {
            const V t0 = v3<S>(0.2, 0, -0.125);
            const V t1 = v3<S>(0.2, 0.25, 0.125);
            const V t2 = v3<S>(0.4, 0, 0);
            if (!hits(a, a, r1, r0, t0, t1, t2)) return fail(__LINE__);

            const V g0 = v3<S>(0.5, 0, -0.125);
            const V g1 = v3<S>(0.5, 0.25, 0.125);
            const V g2 = v3<S>(0.75, 0, 0);
            if (hits(a, a, r1, r0, g0, g1, g2)) return fail(__LINE__);
        }

        /* The degenerate faces: an edge, then a point — the nearest point of the triangle
        holds for both, and the beam answers them as what they are. */
        {
            const V t0 = v3<S>(0.0625, 0, 0.5);
            const V t1 = v3<S>(0.0625, 0, 0.5);
            const V t2 = v3<S>(0.5, 0, 0.5);
            if (!hits(a, b, r0, r1, t0, t1, t2)) return fail(__LINE__);

            const V g0 = v3<S>(0.3, 0, 0.5);
            const V g1 = v3<S>(0.3, 0, 0.5);
            const V g2 = v3<S>(0.3, 0, 0.5);
            if (hits(a, b, r0, r1, g0, g1, g2)) return fail(__LINE__);

            const V p0 = v3<S>(0.0625, 0, 0.5);
            if (!hits(a, b, r0, r1, p0, p0, p0)) return fail(__LINE__);
        }

        /* The nearest point itself: the front face, a vertex's region, an edge's region. */
        {
            const V tri0 = v3<S>(0, -1, 0);
            const V tri1 = v3<S>(0, 1, 0);
            const V tri2 = v3<S>(0, 0, 1);

            auto offBy = [](const V& one, const V& other)
            {
                V delta = one - other;
                return delta.getSquaredNorm();
            };

            const V front = MathGeometry::closestPointOnTriangle(v3<S>(0.25, 0.1, 0.5), tri0, tri1, tri2);
            if (offBy(front, v3<S>(0, 0.1, 0.5)) > S(1e-12)) return fail(__LINE__);

            const V corner = MathGeometry::closestPointOnTriangle(v3<S>(1, 2, 3), tri0, tri1, tri2);
            if (offBy(corner, v3<S>(0, 0, 1)) > S(1e-12)) return fail(__LINE__);

            const V along = MathGeometry::closestPointOnTriangle(v3<S>(1, 0.6, 0.6), tri0, tri1, tri2);
            if (offBy(along, v3<S>(0, 0.5, 0.5)) > S(1e-12)) return fail(__LINE__);
        }

        return true;
    }

bool MathBeamCapsuleTest(std::string& )
{
    if (!checkBeamCapsule<double>())
        return false;
    if (!checkBeamCapsule<float>())
        return false;

    std::cout << "MathBeamCapsuleTest SUCCESS" << std::endl;

    return true;
}

/* The geometry the shadow discontinuities are traced with: the subsets of a facets description, the roots
of a quadratic, the meeting point of two lines, the frame of a plane and the conic swept by the transversals
of three lines. */
bool MathShadowGeometryTest(std::string&)
{
#define SHADOW_CHECK(condition) if (!(condition)) { std::cout << "Error in line " << __LINE__ << ": " #condition << std::endl; return false; }

	{//MathCombinatorial::forEachSubset
		std::vector<std::vector<size_t>> subsets;
		MathCombinatorial::forEachSubset(std::vector<size_t>{ 2, 5, 7, 9 }, 3, [&](const std::vector<size_t>& s) { subsets.push_back(s); });
		SHADOW_CHECK(subsets.size() == 4);
		SHADOW_CHECK((subsets[0] == std::vector<size_t>{ 2, 5, 7 }));
		SHADOW_CHECK((subsets[1] == std::vector<size_t>{ 2, 5, 9 }));
		SHADOW_CHECK((subsets[2] == std::vector<size_t>{ 2, 7, 9 }));
		SHADOW_CHECK((subsets[3] == std::vector<size_t>{ 5, 7, 9 }));

		size_t count = 0;
		MathCombinatorial::forEachSubset(std::vector<size_t>{ 1, 2, 3, 4, 5, 6 }, 3, [&](const std::vector<size_t>&) { count++; });
		SHADOW_CHECK(count == 20);

		count = 0;
		MathCombinatorial::forEachSubset(std::vector<size_t>{ 1, 2 }, 3, [&](const std::vector<size_t>&) { count++; });
		SHADOW_CHECK(count == 0);
		MathCombinatorial::forEachSubset(std::vector<size_t>{ 1, 2 }, 0, [&](const std::vector<size_t>&) { count++; });
		SHADOW_CHECK(count == 0);
	}

	{//MathArithmetic::solveHalfQuadratic: a t^2 + 2 b t + c
		const double tol = 1e-12;
		std::vector<double> roots;
		// (t - 1)(t - 3) = t^2 - 4 t + 3
		SHADOW_CHECK(MathArithmetic<double>::solveHalfQuadratic(1.0, -2.0, 3.0, tol, roots) == 2);
		std::sort(roots.begin(), roots.end());
		SHADOW_CHECK(std::fabs(roots[0] - 1.0) < 1e-12 && std::fabs(roots[1] - 3.0) < 1e-12);

		roots.clear();
		// (t - 2)^2 = t^2 - 4 t + 4: one double root
		SHADOW_CHECK(MathArithmetic<double>::solveHalfQuadratic(1.0, -2.0, 4.0, tol, roots) == 1);
		SHADOW_CHECK(std::fabs(roots[0] - 2.0) < 1e-12);

		roots.clear();
		// t^2 + 1: no real root
		SHADOW_CHECK(MathArithmetic<double>::solveHalfQuadratic(1.0, 0.0, 1.0, tol, roots) == 0);

		roots.clear();
		// linear: 2 * 3 t - 6 = 0
		SHADOW_CHECK(MathArithmetic<double>::solveHalfQuadratic(0.0, 3.0, -6.0, tol, roots) == 1);
		SHADOW_CHECK(std::fabs(roots[0] - 1.0) < 1e-12);

		double t = MathArithmetic<double>::refineHalfQuadraticRoot(1.0, -2.0, 3.0, 1.1, 1e-15);
		SHADOW_CHECK(std::fabs(t - 1.0) < 1e-12);
	}

	{//MathPolynomial::findRoots
		// (t - 0.1)(t - 0.5)(t - 0.9) = t^3 - 1.5 t^2 + 0.59 t - 0.045
		std::vector<double> roots;
		MathPolynomiald({ -0.045, 0.59, -1.5, 1.0 }).findRoots(0.0, 1.0, 1e-12, roots);
		SHADOW_CHECK(roots.size() == 3);
		SHADOW_CHECK(std::fabs(roots[0] - 0.1) < 1e-12 && std::fabs(roots[1] - 0.5) < 1e-12 && std::fabs(roots[2] - 0.9) < 1e-12);

		// only the roots in the interval
		roots.clear();
		MathPolynomiald({ -0.045, 0.59, -1.5, 1.0 }).findRoots(0.3, 1.0, 1e-12, roots);
		SHADOW_CHECK(roots.size() == 2);

		// a double root at a critical point: (t - 0.5)^2 (t + 1)
		roots.clear();
		MathPolynomiald({ 0.25, -0.75, 0.0, 1.0 }).findRoots(0.0, 1.0, 1e-12, roots);
		SHADOW_CHECK(roots.size() == 1 && std::fabs(roots[0] - 0.5) < 1e-6);

		// a null leading coefficient lowers the degree; a null polynomial has no root
		roots.clear();
		MathPolynomiald({ -0.5, 1.0, 0.0, 0.0 }).findRoots(0.0, 1.0, 1e-12, roots);
		SHADOW_CHECK(roots.size() == 1 && std::fabs(roots[0] - 0.5) < 1e-12);
		roots.clear();
		MathPolynomiald({ 0.0, 0.0 }).findRoots(0.0, 1.0, 1e-12, roots);
		SHADOW_CHECK(roots.empty());
	}

	{//MathPlucker6::meets / getMeetingPoint
		const MathVector3d p(1.0, 2.0, 3.0);
		MathPlucker6d l1(p, p + MathVector3d(1.0, 0.5, 0.0));
		MathPlucker6d l2(p + MathVector3d(0.0, -1.0, 2.0), p);
		MathPlucker6d skew(MathVector3d(0.0, 0.0, 10.0), MathVector3d(0.0, 1.0, 10.0));
		MathPlucker6d parallel(p + MathVector3d(0.0, 0.0, 1.0), p + MathVector3d(1.0, 0.5, 1.0));

		SHADOW_CHECK(l1.meets(l2, 1e-9));
		SHADOW_CHECK(!l1.meets(skew, 1e-9));
		SHADOW_CHECK(!l1.meets(parallel, 1e-9));

		MathVector3d meeting;
		SHADOW_CHECK(l1.getMeetingPoint(l2, meeting, 1e-12));
		SHADOW_CHECK((meeting - p).getNorm() < 1e-9);
		SHADOW_CHECK(!l1.getMeetingPoint(parallel, meeting, 1e-12));

		MathPlucker6d l3 = MathPlucker6d::fromPointAndDirection(p, MathVector3d(1.0, 0.5, 0.0));
		SHADOW_CHECK((l3.getDirection() - l1.getDirection()).getNorm() < 1e-12);
		SHADOW_CHECK((l3.getLocation() - l1.getLocation()).getNorm() < 1e-12);
	}

	{//MathFrame3
		MathFrame3d frame(MathVector3d(1.0, 1.0, 2.0), MathVector3d(0.0, 0.0, 3.0), MathVector3d(1.0, 1.0, 1.0));
		SHADOW_CHECK(std::fabs(frame.getAxisU().getNorm() - 1.0) < 1e-12);
		SHADOW_CHECK(std::fabs(frame.getAxisU().dot(frame.getNormal())) < 1e-12);
		SHADOW_CHECK((MathVector3d::cross(frame.getAxisU(), frame.getAxisV()) - frame.getNormal()).getNorm() < 1e-12);

		const MathVector2d uv(0.3, -1.7);
		const MathVector3d world = frame.toWorld(uv);
		SHADOW_CHECK(std::fabs(frame.getSignedDistance(world)) < 1e-12);
		SHADOW_CHECK((frame.toLocal(world) - uv).getNorm() < 1e-12);

		const MathPlane3d plane = frame.getPlane();
		SHADOW_CHECK(std::fabs(plane.mNormal.dot(world) + plane.d) < 1e-12);
	}

	{//MathConic2::fromTransversalsOfThreeLines
		// three skew lines, and their transversals: every transversal crosses the plane z = 0 on the conic
		const MathVector3d a[3] = { MathVector3d(0.0, 0.0, 1.0), MathVector3d(0.0, 0.0, 2.0), MathVector3d(0.3, 0.2, 3.0) };
		const MathVector3d d[3] = { MathVector3d(1.0, 0.0, 0.0), MathVector3d(0.0, 1.0, 0.0), MathVector3d(1.0, 1.0, 0.2) };
		MathPlucker6d lines[3];
		for (int i = 0; i < 3; i++)
		{
			lines[i] = MathPlucker6d::fromPointAndDirection(a[i], d[i]);
		}
		MathFrame3d frame(MathVector3d(0.1, -0.2, 0.0), MathVector3d(0.0, 0.0, 1.0), MathVector3d(1.0, 0.2, 0.0));
		MathConic2d conic = MathConic2d::fromTransversalsOfThreeLines(lines, frame);
		SHADOW_CHECK(!conic.isLinear(1e-9));

		size_t tested = 0;
		for (double t = -2.0; t <= 2.0; t += 0.25)
		{
			// the transversal through a(t) on line 0: it meets line 2 where line 2 crosses the plane (a(t), line 1)
			const MathVector3d at = a[0] + d[0] * t;
			MathVector3d normal = MathVector3d::cross(d[1], a[1] - at);
			const double denominator = normal.dot(d[2]);
			if (std::fabs(denominator) < 1e-9)
			{
				continue;
			}
			const double s = normal.dot(at - a[2]) / denominator;
			const MathVector3d ct = a[2] + d[2] * s;
			const MathVector3d direction = ct - at;
			if (std::fabs(direction.z) < 1e-9)
			{
				continue;
			}
			const MathVector3d onPlane = at - direction * (at.z / direction.z);

			const MathPlucker6d transversal(at, ct);
			SHADOW_CHECK(std::fabs(transversal.dot(lines[1])) < 1e-9);
			SHADOW_CHECK(conic.getAlgebraicDistance(frame.toLocal(onPlane)) < 1e-9);
			tested++;
		}
		SHADOW_CHECK(tested > 8);

		// a point away from the regulus is not on the conic
		SHADOW_CHECK(conic.getAlgebraicDistance(MathVector2d(10.0, -7.0)) > 1e-3);
	}
#undef SHADOW_CHECK
	return true;
}

/* The incidences of the vertices with the hyperplanes, decided from the facet descriptions, and the merge of the
hyperplanes carrying the same line. */
bool MathIncidenceTest(std::string&)
{
#define INCIDENCE_CHECK(condition) if (!(condition)) { std::cout << "Error in line " << __LINE__ << ": " #condition << std::endl; return false; }
	const double tolerance = MathArithmetic<double>::Tolerance();

	{//MathCombinatorial::initFacets of a merge: the union of the two descriptions
		std::vector<size_t> result;
		MathCombinatorial::initFacets(std::vector<size_t>{ 1, 2, 5 }, std::vector<size_t>{ 2, 3, 5 }, result);
		INCIDENCE_CHECK((result == std::vector<size_t>{ 1, 2, 3, 5 }));

		result.clear();
		MathCombinatorial::initFacets(std::vector<size_t>{ 4, 7 }, std::vector<size_t>{ 4, 7 }, result);
		INCIDENCE_CHECK((result == std::vector<size_t>{ 4, 7 }));
	}

	{//MathPredicates::getProportionality
		const MathPlucker6d a(MathVector3d(0.1, 0.2, 0.3), MathVector3d(1.0, -0.5, 2.0));
		const MathPlucker6d other(MathVector3d(0.1, 0.2, 0.3), MathVector3d(1.0, -0.5, 2.1));
		INCIDENCE_CHECK(MathPredicates::getProportionality(a, a * 3.0, tolerance) == 1);
		INCIDENCE_CHECK(MathPredicates::getProportionality(a, a * -2.0, tolerance) == -1);
		INCIDENCE_CHECK(MathPredicates::getProportionality(a, other, tolerance) == 0);

		// a rounding far below the tolerance is the same plane, a difference above it is not
		MathPlucker6d nearly = a * 2.0;
		nearly += MathPlucker6d(1e-15, -1e-15, 0.0, 0.0, 1e-15, 0.0);
		INCIDENCE_CHECK(MathPredicates::getProportionality(a, nearly, tolerance) == 1);
		MathPlucker6d apart = a;
		apart += MathPlucker6d(0.0, 1e-6, 0.0, 0.0, 0.0, 0.0);
		INCIDENCE_CHECK(MathPredicates::getProportionality(a, apart, tolerance) == 0);
	}

	{//PluckerPolyhedron::addHyperplane: the same line, either way, is merged; another line is not
		PluckerPolyhedron<MathPlucker6d> polyhedron;
		const MathPlucker6d a(MathVector3d(0.0, 0.0, 0.0), MathVector3d(1.0, 0.0, 0.0));
		const MathPlucker6d collinear(MathVector3d(3.0, 0.0, 0.0), MathVector3d(2.0, 0.0, 0.0)); // same line, the other way
		const MathPlucker6d other(MathVector3d(0.0, 1.0, 0.0), MathVector3d(1.0, 1.0, 0.0));

		int orientation = 0;
		bool merged = true;
		const size_t first = polyhedron.addHyperplane(a, false, tolerance, orientation, merged);
		INCIDENCE_CHECK(!merged && orientation == 1);

		const size_t second = polyhedron.addHyperplane(collinear, false, tolerance, orientation, merged);
		INCIDENCE_CHECK(merged && second == first && orientation == -1);

		const size_t third = polyhedron.addHyperplane(other, false, tolerance, orientation, merged);
		INCIDENCE_CHECK(!merged && third != first && orientation == 1);

		// a plain add is never a candidate for a merge
		const size_t plain = polyhedron.add(other * 5.0, ON_BOUNDARY, false, tolerance);
		INCIDENCE_CHECK(polyhedron.addHyperplane(other, false, tolerance, orientation, merged) == third && plain != third);
	}

	{//SilhouetteEdge: a side of an edge merged the other way is the opposite side of its hyperplane
		SilhouetteEdge edge(nullptr, 0, 0);
		INCIDENCE_CHECK(edge.getEdgeSide(ON_NEGATIVE_SIDE) == ON_NEGATIVE_SIDE);
		edge.mHyperPlaneSign = -1;
		INCIDENCE_CHECK(edge.getEdgeSide(ON_NEGATIVE_SIDE) == ON_POSITIVE_SIDE);
		INCIDENCE_CHECK(edge.getEdgeSide(ON_POSITIVE_SIDE) == ON_NEGATIVE_SIDE);
		INCIDENCE_CHECK(edge.getEdgeSide(ON_BOUNDARY) == ON_BOUNDARY);
		INCIDENCE_CHECK(edge.getEdgeSide(ON_BOTH_SIDES) == ON_BOTH_SIDES);
	}

	{//the splitter reads the incidence of a vertex with one of its facets from its description
		// two triangles facing each other: the initial cell, whose vertices lie on the edges of the two sources
		const float light[] = { -1.0f, 0.0f, 0.3f, -1.0f, 0.26f, -0.15f, -1.0f, -0.26f, -0.15f };
		const float receiver[] = { 1.0f, 0.0f, 0.3f, 1.0f, -0.26f, -0.15f, 1.0f, 0.26f, -0.15f };
		GeometryConvexPolygon a(light, 3);
		GeometryConvexPolygon b(receiver, 3);

		PluckerPolyhedron<MathPlucker6d> polyhedron;
		PluckerPolytopeBuilder<MathPlucker6d, double> builder(MathPluckerNormalization<MathPlucker6d, double>(), tolerance);
		PluckerPolytope<MathPlucker6d>* polytope = builder.build(a, b, &polyhedron);
		INCIDENCE_CHECK(polytope->getEdgeCount() > 0);

		// the hyperplane of the first edge of the light, a facet of the cell, read with an error a thousand times the
		// tolerance: the vertices lying on it must still be found on it, and the cell on one side, unsplit
		/* An error that separates two vertices of the facet: with u1, u2 the two vertices on the unit sphere, the swap of
		u1 - u2 gives a Klein product with a vertex v equal to (u1 - u2) . v, positive on the first and negative on the
		second (Cauchy-Schwarz). */
		std::vector<size_t> incident;
		for (size_t vertex : polytope->getVertices())
		{
			if (MathCombinatorial::hasFacet(polyhedron.getFacetsDescription(vertex), 0))
			{
				incident.push_back(vertex);
			}
		}
		INCIDENCE_CHECK(incident.size() >= 2);
		const MathPlucker6d u1 = polyhedron.get(incident[0]).getNormalized();
		const MathPlucker6d u2 = polyhedron.get(incident[1]).getNormalized();
		MathPlucker6d plane = polyhedron.get(0);
		plane += (u1 - u2).quadricGradient() * 1e-8;

		// the error must put the incident vertices on both sides numerically, or the check below proves nothing
		bool hasNegative = false;
		bool hasPositive = false;
		for (size_t vertex : polytope->getVertices())
		{
			if (MathCombinatorial::hasFacet(polyhedron.getFacetsDescription(vertex), 0))
			{
				const GeometryPositionType numeric = MathPredicates::getVertexPlaneRelativePosition(plane, polyhedron.get(vertex), tolerance);
				hasNegative |= numeric == ON_NEGATIVE_SIDE;
				hasPositive |= numeric == ON_POSITIVE_SIDE;
			}
		}
		INCIDENCE_CHECK(hasNegative && hasPositive);

		const size_t pointsBefore = polyhedron.getLinesCount();
		PluckerPolytope<MathPlucker6d>* left = new PluckerPolytope<MathPlucker6d>();
		PluckerPolytope<MathPlucker6d>* right = new PluckerPolytope<MathPlucker6d>();
		HelperStatisticCollector statistic;
		const GeometryPositionType result = PluckerPolytopeSplitter<MathPlucker6d, double>::split(&polyhedron, plane, polytope, left, right, 0, MathPluckerNormalization<MathPlucker6d, double>(), tolerance, &statistic);

		INCIDENCE_CHECK(result != ON_BOTH_SIDES);
		INCIDENCE_CHECK(polyhedron.getLinesCount() == pointsBefore);
		INCIDENCE_CHECK(statistic.get(INCIDENCE_COMBINATORIAL_COUNT) > 0);

		delete left;
		delete right;
		delete polytope;
	}
#undef INCIDENCE_CHECK
	return true;
}

/* The exact predicates on the points of an integer grid: the 128 bit products, the snapping, the orientation of four
points, and its agreement with the Klein product of the Plucker lines through them. */
bool MathIntegerPredicateTest(std::string&)
{
#define INTEGER_CHECK(condition) if (!(condition)) { std::cout << "Error in line " << __LINE__ << ": " #condition << std::endl; return false; }

	{//MathInt128
		const int64_t big = int64_t(1) << 62;
		INTEGER_CHECK(MathInt128::multiply(big, big).sign() == 1);
		INTEGER_CHECK(MathInt128::multiply(big, -big).sign() == -1);
		INTEGER_CHECK((MathInt128::multiply(big, big) - MathInt128::multiply(big, big)).sign() == 0);
		// 2^124 - (2^124 - 1): the difference is in the lowest bit
		const MathInt128 a = MathInt128::multiply(big, big);
		const MathInt128 b = MathInt128::multiply(big, big) - MathInt128(1);
		INTEGER_CHECK((a - b).sign() == 1);
		INTEGER_CHECK((b - a).sign() == -1);
		INTEGER_CHECK((MathInt128::multiply(-3, 5) + MathInt128(15)).sign() == 0);
	}

	{//MathIntegerGrid
		MathIntegerGrid grid(3.0);
		INTEGER_CHECK(grid.isValid());
		INTEGER_CHECK(grid.getExtent() >= 3.0 && grid.getExtent() < 6.0 + 1e-12);
		const MathVector3f p(0.123456789f, -2.5f, 1e-9f);
		const MathVector3f q = grid.snap(p);
		INTEGER_CHECK(grid.isOnGrid(q));
		// the floats of ordinary magnitude are already on the grid: only the tiny ones move
		INTEGER_CHECK(q.x == p.x && q.y == p.y);
		INTEGER_CHECK(q.z != p.z && std::fabs(q.z - p.z) <= 0.5 * grid.getStep());
		INTEGER_CHECK(std::fabs(q.x - p.x) <= 0.5 * grid.getStep() && std::fabs(q.z - p.z) <= 0.5 * grid.getStep());
		INTEGER_CHECK(!grid.isOnGrid(MathVector3d(0.1 * grid.getStep(), 0.0, 0.0)));
		INTEGER_CHECK(!grid.isOnGrid(MathVector3d(4.0 * grid.getExtent(), 0.0, 0.0)));
		const MathIntegerPoint i = grid.toInteger(q);
		INTEGER_CHECK(double(i.y) * grid.getStep() == double(q.y));
		INTEGER_CHECK(!MathIntegerGrid().isValid());
	}

	{//orient3d, kleinSign and areCollinear
		MathIntegerGrid grid(1.0);
		srand(7);
		auto randomPoint = [&]()
		{
			return grid.snap(MathVector3f(float(rand()) / RAND_MAX * 2.0f - 1.0f, float(rand()) / RAND_MAX * 2.0f - 1.0f, float(rand()) / RAND_MAX * 2.0f - 1.0f));
		};
		size_t compared = 0;
		for (int k = 0; k < 2000; k++)
		{
			const MathVector3f a = randomPoint(), b = randomPoint(), c = randomPoint(), d = randomPoint();
			const MathIntegerPoint ia = grid.toInteger(a), ib = grid.toInteger(b), ic = grid.toInteger(c), id = grid.toInteger(d);

			// the Klein product of the two lines, in double: its sign must be the exact one whenever it is clearly away from zero
			const double klein = MathPlucker6d(convert<MathVector3d>(a), convert<MathVector3d>(b)).dot(MathPlucker6d(convert<MathVector3d>(c), convert<MathVector3d>(d)));
			if (std::fabs(klein) > 1e-9)
			{
				INTEGER_CHECK(MathIntegerPredicates::kleinSign(ia, ib, ic, id) == (klein > 0 ? 1 : -1));
				compared++;
			}

			// a fourth point built in the plane of the three first, exactly: the orientation is exactly zero
			const MathIntegerPoint planar{ ia.x + 2 * (ib.x - ia.x) - 3 * (ic.x - ia.x), ia.y + 2 * (ib.y - ia.y) - 3 * (ic.y - ia.y), ia.z + 2 * (ib.z - ia.z) - 3 * (ic.z - ia.z) };
			INTEGER_CHECK(MathIntegerPredicates::orient3d(ia, ib, ic, planar) == 0);

			// a third point on the line of the two first
			const MathIntegerPoint aligned{ ia.x + 5 * (ib.x - ia.x), ia.y + 5 * (ib.y - ia.y), ia.z + 5 * (ib.z - ia.z) };
			INTEGER_CHECK(MathIntegerPredicates::areCollinear(ia, ib, aligned));
			INTEGER_CHECK(MathIntegerPredicates::orient3d(ia, ib, aligned, id) == 0);
		}
		INTEGER_CHECK(compared > 1900);

		// the largest grid coordinates: the products reach the top of the 128 bits without overflow
		const int64_t m = (int64_t(1) << (MathIntegerGrid::cBits - 1));
		const MathIntegerPoint o{ -m, -m, -m }, x{ m, -m, -m }, y{ -m, m, -m }, z{ -m, -m, m };
		INTEGER_CHECK(MathIntegerPredicates::orient3d(o, x, y, z) == 1);
		INTEGER_CHECK(MathIntegerPredicates::orient3d(o, y, x, z) == -1);
	}
#undef INTEGER_CHECK
	return true;
}

/* The three normalizations of the Plucker points: none, the unit hypersphere, and the chart of the source and the
receiver. A normalization multiplies a point by a positive factor of its own, so no Klein sign may change; the chart
puts the lines from the source to the receiver, and their positive combinations, at l = 1, and scales the hyperplanes
by a power of two. */
bool MathNormalizationTest(std::string&)
{
#define NORMALIZATION_CHECK(condition) if (!(condition)) { std::cout << "Error in line " << __LINE__ << ": " #condition << std::endl; return false; }
	typedef MathPluckerNormalization<MathPlucker6d, double> Normalization;

	// a light and a receiver facing each other, slightly tilted: the receiver on the positive side of the light plane
	const float light[] = { -1.0f, 0.0f, 0.3f, -1.1f, 0.26f, -0.15f, -0.9f, -0.26f, -0.15f };
	const float receiver[] = { 1.0f, 0.0f, 0.3f, 1.2f, -0.26f, -0.15f, 0.9f, 0.26f, -0.15f };
	GeometryConvexPolygon a(light, 3);
	GeometryConvexPolygon b(receiver, 3);
	if (a.getPlane().dot(b.getVertex(0)) < 0.0)
	{
		// the windings decide the planes: turn them so that each polygon faces the other, as the query does
		a = GeometryConvexPolygon(std::vector<MathVector3d>{ a.getVertex(0), a.getVertex(2), a.getVertex(1) });
	}
	if (b.getPlane().dot(a.getVertex(0)) < 0.0)
	{
		b = GeometryConvexPolygon(std::vector<MathVector3d>{ b.getVertex(0), b.getVertex(2), b.getVertex(1) });
	}

	double margin = 0.0;
	const MathVector3d chartNormal = Normalization::computeChartNormal(a, b, margin);
	NORMALIZATION_CHECK(margin > 0.5);
	NORMALIZATION_CHECK(std::fabs(chartNormal.getNorm() - 1.0) < 1e-12);

	const Normalization none(NO_NORMALIZATION, chartNormal);
	const Normalization sphere(HYPERSPHERE_NORMALIZATION, chartNormal);
	const Normalization chart(CHART_NORMALIZATION, chartNormal);
	const Normalization* modes[3] = { &none, &sphere, &chart };

	{//the stabbing lines and their positive combinations land in the chart
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				const MathPlucker6d line(a.getVertex(i), b.getVertex(j));
				NORMALIZATION_CHECK(chart.getChartValue(line) > 0.0);
				const MathPlucker6d inChart = chart.normalizePoint(line);
				NORMALIZATION_CHECK(std::fabs(chart.getChartValue(inChart) - 1.0) < 1e-12);

				const MathPlucker6d other = chart.normalizePoint(MathPlucker6d(a.getVertex((i + 1) % 3), b.getVertex((j + 2) % 3)));
				const MathPlucker6d combination = inChart * 0.3 + other * 0.7;
				NORMALIZATION_CHECK(std::fabs(chart.getChartValue(combination) - 1.0) < 1e-12);
			}
		}
	}

	{//a reversed line keeps its orientation: divided by |l|, not by l
		const MathPlucker6d forward(a.getVertex(0), b.getVertex(1));
		const MathPlucker6d backward(b.getVertex(1), a.getVertex(0));
		const MathPlucker6d edge(MathVector3d(0.0, -0.3, 0.1), MathVector3d(0.1, 0.4, -0.2));
		const MathPlucker6d forwardInChart = chart.normalizePoint(forward);
		const MathPlucker6d backwardInChart = chart.normalizePoint(backward);
		NORMALIZATION_CHECK(chart.getChartValue(backwardInChart) < 0.0);
		NORMALIZATION_CHECK((forwardInChart.dot(edge) > 0.0) == (forward.dot(edge) > 0.0));
		NORMALIZATION_CHECK((backwardInChart.dot(edge) > 0.0) == (backward.dot(edge) > 0.0));
	}

	{//the hyperplanes are scaled by a power of two: exactly, and their largest coordinate in [0.5, 1)
		const MathPlucker6d edge(MathVector3d(0.123, -4.56, 7.89), MathVector3d(-3.21, 6.54, -0.987));
		const MathPlucker6d scaled = chart.normalizeHyperplane(edge);
		const double c[6] = { scaled.getDirection().x, scaled.getDirection().y, scaled.getDirection().z, scaled.getLocation().x, scaled.getLocation().y, scaled.getLocation().z };
		const double o[6] = { edge.getDirection().x, edge.getDirection().y, edge.getDirection().z, edge.getLocation().x, edge.getLocation().y, edge.getLocation().z };
		double largest = 0.0;
		for (double v : c) largest = std::fmax(largest, std::fabs(v));
		NORMALIZATION_CHECK(largest >= 0.5 && largest < 1.0);
		const double ratio = c[0] / o[0];
		int exponent = 0;
		NORMALIZATION_CHECK(std::frexp(ratio, &exponent) == 0.5);
		for (int k = 0; k < 6; k++)
		{
			NORMALIZATION_CHECK(c[k] == o[k] * ratio); // exact: no rounding
		}
		// an edge parallel to the chart plane, which has no chart value, is a hyperplane like any other
		const MathPlucker6d parallel = MathPlucker6d::fromPointAndDirection(MathVector3d(0.0, 0.2, 0.1), MathVector3d::cross(chartNormal, MathVector3d(0.0, 0.0, 1.0)));
		NORMALIZATION_CHECK(std::fabs(chart.getChartValue(parallel)) < 1e-12);
		NORMALIZATION_CHECK(chart.normalizeHyperplane(parallel).getDirection().getNorm() > 0.0);
	}

	{//no Klein sign depends on the mode
		srand(11);
		auto randomPoint = []() { return MathVector3d(rand() / (double)RAND_MAX - 0.5, rand() / (double)RAND_MAX - 0.5, rand() / (double)RAND_MAX - 0.5); };
		for (int k = 0; k < 500; k++)
		{
			const MathPlucker6d line(a.getVertex(k % 3) + randomPoint() * 0.1, b.getVertex((k + 1) % 3) + randomPoint() * 0.1);
			const MathPlucker6d edge(randomPoint(), randomPoint());
			const double reference = line.dot(edge);
			if (std::fabs(reference) < 1e-9)
			{
				continue;
			}
			for (const Normalization* mode : modes)
			{
				const double value = mode->normalizePoint(line).dot(mode->normalizeHyperplane(edge));
				NORMALIZATION_CHECK((value > 0.0) == (reference > 0.0));
			}
		}
	}
#undef NORMALIZATION_CHECK
	return true;
}

/* Simulation of Simplicity in GMP rational: every sign of the perturbed scene is exact and never zero unless the
incidence is structural, the constructions carry their terms of epsilon, and a split by a hyperplane exactly through
vertices of the cell stays consistent. */
bool MathSimulationOfSimplicityTest(std::string&)
{
#ifdef ENABLE_GMP
#define SOS_CHECK(condition) if (!(condition)) { std::cout << "Error in line " << __LINE__ << ": " #condition << std::endl; return false; }
	typedef MathGmpRational Q;
	typedef MathPlucker6<Q> L;
	typedef MathSoS<L, Q> SoS;
	const Q tolerance = MathArithmetic<Q>::Tolerance();
	auto point = [](double x, double y, double z) { return MathVector3d(x, y, z); };

	{//the perturbation directions: deterministic, never null, different from one index to the next
		for (size_t i = 0; i < 200; i++)
		{
			SOS_CHECK(!SoS::isZero(SoS::perturbationOf(i)));
			SOS_CHECK(SoS::isZero(SoS::perturbationOf(i) - SoS::perturbationOf(i)));
		}
		size_t different = 0;
		for (size_t i = 0; i + 1 < 200; i++)
		{
			different += SoS::isZero(SoS::perturbationOf(i) - SoS::perturbationOf(i + 1)) ? 0 : 1;
		}
		SOS_CHECK(different == 199);
	}

	// a vertex: the line through a and b; a hyperplane through a meets it, exactly
	const MathVector3d a = point(-1.0, 0.25, 0.5), b = point(1.0, -0.5, 0.125), c = point(0.5, 2.0, -1.0);
	const L vertexLine(a, b);
	const L plane(a, c);
	const L u = SoS::perturbationOf(42);
	SOS_CHECK(getExactSign(plane.dot(vertexLine)) == 0);

	{//an accidental zero is broken at order 1, the same way every time; a structural one is a residual
		size_t order = 0;
		const int s1 = SoS::sign(plane, &u, MathSymbolicPoint<L>(vertexLine), order);
		SOS_CHECK(s1 != 0 && order == 1);
		const int s2 = SoS::sign(plane, &u, MathSymbolicPoint<L>(vertexLine), order);
		SOS_CHECK(s2 == s1);
		// the same plane unperturbed: nothing breaks the tie
		SOS_CHECK(SoS::sign(plane, nullptr, MathSymbolicPoint<L>(vertexLine), order) == 0 && order == SoS::cOrder + 1);
	}

	{//a crossing at an end whose value vanishes at order 0: the leading term is that end, the point another one
		size_t order = 0;
		const int sEnd = SoS::sign(plane, &u, MathSymbolicPoint<L>(vertexLine), order);
		// the other end: a line strictly on the other side of the plane, at order 0
		const L candidates[4] = { L(point(-1.0, 0.75, -0.5), point(1.0, 0.5, 0.5)), L(point(1.0, 0.5, 0.5), point(-1.0, 0.75, -0.5)),
		                          L(point(-1.0, -0.75, 0.5), point(1.0, 0.25, -0.5)), L(point(1.0, 0.25, -0.5), point(-1.0, -0.75, 0.5)) };
		L other = candidates[0];
		for (const L& candidate : candidates)
		{
			if (getExactSign(plane.dot(candidate)) == -sEnd)
			{
				other = candidate;
				break;
			}
		}
		const int sOther = SoS::sign(plane, &u, MathSymbolicPoint<L>(other), order);
		SOS_CHECK(sEnd != 0 && sOther != 0);
		SOS_CHECK(sEnd != sOther); // the edge crosses the perturbed plane: the case under test
		{
			const MathSymbolicPoint<L> negative(sEnd < 0 ? vertexLine : other);
			const MathSymbolicPoint<L> positive(sEnd < 0 ? other : vertexLine);
			const MathSymbolicPoint<L> crossing = SoS::crossing(plane, &u, negative, positive);
			// the leading term is a positive multiple of the end on the plane
			SOS_CHECK(MathPredicates::getProportionality(crossing.leading(), vertexLine, Q(0)) == 1);
			SOS_CHECK(!SoS::areEqual(crossing, MathSymbolicPoint<L>(vertexLine)));
			SOS_CHECK(crossing.size() > 1);
			// it lies on the perturbed plane: every coefficient of the plane at the crossing vanishes
			SOS_CHECK(SoS::sign(plane, &u, crossing, order) == 0);
		}
	}

	{//the quadric: a line is real by construction (a residual), an interpolated point of two skew lines is not
		size_t order = 0;
		SOS_CHECK(SoS::quadricSign(MathSymbolicPoint<L>(vertexLine), order) == 0);
		const L skew(point(0.0, 1.0, 0.0), point(0.0, 1.0, 1.0));
		SOS_CHECK(getExactSign(vertexLine.dot(skew)) != 0);
		SOS_CHECK(SoS::quadricSign(MathSymbolicPoint<L>(vertexLine + skew), order) != 0);
	}

	{//the exact quadric crossing of an edge agrees with the double one on generic edges
		srand(5);
		auto r = []() { return rand() / (double)RAND_MAX - 0.5; };
		size_t compared = 0;
		for (int k = 0; k < 400; k++)
		{
			const MathVector3d p0 = point(r(), r(), r()), p1 = point(r(), r(), r()), p2 = point(r(), r(), r()), p3 = point(r(), r(), r());
			const L v1 = L(p0, p1) + L(p2, p3) * Q(r());
			const L v2 = L(p1, p2) + L(p3, p0) * Q(r());
			size_t order = 0;
			const int q1 = SoS::quadricSign(MathSymbolicPoint<L>(v1), order);
			const int q2 = SoS::quadricSign(MathSymbolicPoint<L>(v2), order);
			const GeometryPositionType e1 = q1 > 0 ? ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE;
			const GeometryPositionType e2 = q2 > 0 ? ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE;
			bool tangent = false;
			const bool exact = SoS::edgeCrossesQuadric(v1, v2, e1, e2, tangent);
			const MathPlucker6d d1 = convert<MathPlucker6d>(v1), d2 = convert<MathPlucker6d>(v2);
			const bool approximate = MathGeometry::hasPluckerEdgeWithQuadricIntersection(d1, d2, e1, e2, 1e-12);
			SOS_CHECK(!tangent);
			SOS_CHECK(exact == approximate);
			compared++;
		}
		SOS_CHECK(compared == 400);
	}

	{//a split of the initial cell by a perturbed hyperplane exactly through one of its vertices
		const float light[] = { -1.0f, 0.0f, 0.25f, -1.0f, 0.25f, -0.125f, -1.0f, -0.25f, -0.125f };
		const float receiver[] = { 1.0f, 0.0f, 0.25f, 1.0f, -0.25f, -0.125f, 1.0f, 0.25f, -0.125f };
		GeometryConvexPolygon lightPolygon(light, 3);
		GeometryConvexPolygon receiverPolygon(receiver, 3);

		PluckerPolyhedron<L> polyhedron;
		PluckerPolytopeBuilder<L, Q> builder(MathPluckerNormalization<L, Q>(), tolerance);
		PluckerPolytope<L>* polytope = builder.build(lightPolygon, receiverPolygon, &polyhedron);
		polyhedron.setSymbolic(true);

		// a silhouette line through the first light vertex and the interior of the receiver: it meets every VV line of that vertex
		const MathVector3d through = convert<MathVector3d>(lightPolygon.getVertex(0));
		const L silhouette(through, point(0.0, 0.0, 0.125)); // it reaches the receiver at (1, 0, 0), inside
		int orientation = 0;
		bool merged = false;
		const size_t index = polyhedron.addHyperplane(silhouette, false, tolerance, orientation, merged);
		polyhedron.setPerturbation(index, SoS::perturbationOf(index));

		HelperStatisticCollector statistic;
		PluckerPolytope<L>* left = new PluckerPolytope<L>();
		PluckerPolytope<L>* right = new PluckerPolytope<L>();
		const GeometryPositionType result = PluckerPolytopeSplitter<L, Q>::split(&polyhedron, silhouette, polytope, left, right, index, MathPluckerNormalization<L, Q>(), tolerance, &statistic);

		SOS_CHECK(statistic.get(SOS_TIE_BREAK_COUNT) > 0);     // the vertices of the light vertex were on the plane, exactly
		SOS_CHECK(statistic.get(SOS_RESIDUAL_COUNT) == 0);     // and the perturbation decided them all
		SOS_CHECK(statistic.get(INCIDENCE_NUMERIC_BOUNDARY_COUNT) == 0);
		SOS_CHECK(result == ON_BOTH_SIDES); // the plane goes through the cell: the split under test
		{
			// each child on its side of the perturbed plane, every vertex strictly or by its facets
			for (PluckerPolytope<L>* child : { left, right })
			{
				const GeometryPositionType side = MathPredicates::getRelativePosition(child, &polyhedron, silhouette, tolerance, index);
				SOS_CHECK(side == (child == left ? ON_NEGATIVE_SIDE : ON_POSITIVE_SIDE));
			}
		}
		delete left;
		delete right;
		delete polytope;
	}
#undef SOS_CHECK
#endif
	return true;
}

/* The lazy exact rational (MathLazyRational): the values of GMP rational, decided by intervals first.
- random expressions, built in parallel in GMP rational: the double and MPFI intervals contain the exact value, the
  exact value is the one of GMP rational, and the comparisons agree with it;
- the exact zeros (cancellations) are decided by the exact stage;
- a division by an interval holding zero, and an overflow of the double stage, are decided by the later stages;
- a deep chain gets its exact value, which prunes its expression. */
bool MathLazyRationalTest(std::string&)
{
#if defined(ENABLE_GMP) && defined(ENABLE_MPFI)
#define LAZY_CHECK(condition) if (!(condition)) { std::cout << "Error in line " << __LINE__ << ": " #condition << std::endl; return false; }
	typedef MathLazyRational L;
	typedef MathGmpRational Q;

	// -1, 0 or +1 as the rational q is below, equal to or above the double d (an infinity bounds everything)
	auto compareToDouble = [](const Q& q, double d)
	{
		if (std::isinf(d)) return d > 0 ? -1 : 1;
		mpq_t t;
		mpq_init(t);
		mpq_set_d(t, d);
		const int c = mpq_cmp(q.v, t);
		mpq_clear(t);
		return c < 0 ? -1 : (c > 0 ? 1 : 0);
	};
	// the intervals of x contain q, read before the exact value of x is asked (it tightens them)
	auto encloses = [&](const L& x, const Q& q)
	{
		if (compareToDouble(q, x.lower()) < 0 || compareToDouble(q, x.upper()) > 0) return false;
		mpfi_t interval;
		mpfi_init2(interval, L::cIntervalPrecision);
		x.getInterval(interval);
		const bool inside = mpfi_is_inside_q(q.v, interval) != 0;
		mpfi_clear(interval);
		return inside;
	};

	{//random expressions
		srand(11);
		auto leaf = [](L& x, Q& q)
		{
			const int kind = rand() % 3;
			if (kind == 0)
			{
				const double d = (rand() / (double)RAND_MAX - 0.5) * std::pow(10.0, rand() % 7 - 3);
				x = L(d);
				q = Q(d);
			}
			else if (kind == 1)
			{
				q = Q(double(rand() % 19 - 9)) / Q(double(rand() % 7 + 1)); // a rational no double is
				x = L(q);
			}
			else
			{
				const double d = double(rand() % 5 - 2);
				x = L(d);
				q = Q(d);
			}
		};
		size_t checked = 0;
		for (int trial = 0; trial < 300; trial++)
		{
			L x;
			Q q;
			leaf(x, q);
			for (int step = 0; step < 40; step++)
			{
				L y;
				Q r;
				leaf(y, r);
				switch (rand() % 5)
				{
				case 0: x = x + y; q = q + r; break;
				case 1: x = x - y; q = q - r; break;
				case 2: x = x * y; q = q * r; break;
				case 3: if (mpq_sgn(r.v) != 0) { x = x / y; q = q / r; } break;
				default: x = -x; q = -q; break;
				}
				LAZY_CHECK(encloses(x, q));
				// a comparison with a fresh value, against the exact one
				L z;
				Q s;
				leaf(z, s);
				const int expected = mpq_cmp(q.v, s.v);
				const int c = L::compare(x, z);
				LAZY_CHECK(c == (expected < 0 ? -1 : (expected > 0 ? 1 : 0)));
				checked++;
			}
			LAZY_CHECK(mpq_equal(x.exact().v, q.v) != 0);
		}
		LAZY_CHECK(checked == 300 * 40);
	}

	{//exact zeros are decided by the exact stage
		const L a(0.1), b(Q(1.0) / Q(3.0)), c(Q(-2.0) / Q(7.0));
		const size_t before = L::getStatistics().exactDecisions;
		const L z = (a + b) - a - b;
		LAZY_CHECK(z.sign() == 0);
		const L w = a * b * c - c * (b * a);
		LAZY_CHECK(w.sign() == 0);
		LAZY_CHECK(L::getStatistics().exactDecisions >= before + 2);
		// a value next to zero, beyond the double intervals, decided by the MPFI ones
		const L tiny = (a + L(1e-30)) - a;
		LAZY_CHECK(tiny.sign() > 0);
	}

	{//a division by an interval holding zero, an overflow
		const L small = L(1e-300) * L(1e-300); // underflows: its double interval holds zero
		LAZY_CHECK(small.lower() <= 0.0 && small.upper() >= 0.0);
		const L huge = L(1.0) / small;
		LAZY_CHECK(std::isinf(huge.upper()));
		LAZY_CHECK(huge > L(1e300));
		const L big = L(1e300) * L(1e300); // overflows
		LAZY_CHECK(std::isinf(big.upper()));
		LAZY_CHECK(big > L(1e308));
		LAZY_CHECK(-big < L(-1e308));
	}

	{//a deep chain: its exact value is the one of GMP rational, and computing it prunes the expression
		L x(Q(1.0) / Q(3.0));
		Q q = Q(1.0) / Q(3.0);
		const L step(0.75), shift(Q(1.0) / Q(5.0));
		const Q stepQ(0.75), shiftQ = Q(1.0) / Q(5.0);
		for (int i = 0; i < 2000; i++)
		{
			x = x * step + shift;
			q = q * stepQ + shiftQ;
		}
		LAZY_CHECK(x.hasOperands());
		LAZY_CHECK(encloses(x, q));
		LAZY_CHECK(mpq_equal(x.exact().v, q.v) != 0);
		LAZY_CHECK(!x.hasOperands());
		LAZY_CHECK(encloses(x, q));
	}
#undef LAZY_CHECK
#endif
	return true;
}
