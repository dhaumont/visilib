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

#include <cmath>
#include <limits>

#include "math_arithmetic.h"
#include "math_plucker_normalization_type.h"
#include "math_sos.h"
#include "math_vector_3.h"

namespace visilib
{
    /** @brief The normalization of the Plucker points of a query: none, the unit hypersphere, or an affine chart.

    Visilib works in oriented projective space: a Plucker point stands for the ray of its positive multiples, and every
    predicate (the sign of a Klein product, the side of the quadric) and every construction (the vertex where an edge
    crosses a hyperplane, the segment of positive combinations between two vertices) is unchanged when a point is
    multiplied by a positive factor of its own. A normalization is therefore free to choose that factor point by point,
    as long as it is positive: what it changes is the magnitude of the values the tolerance is compared to, and the
    rounding of the constructions.

    - NO_NORMALIZATION keeps the points as they are built.
    - HYPERSPHERE_NORMALIZATION divides each point by its six dimensional norm (a square root).
    - CHART_NORMALIZATION treats the two kinds of points differently:
      - a cell point (a vertex, or a line used as a point: a stabbing line, a representative line) is divided by the
        linear form l(v) = d . n, d its direction part and n the chart normal (see computeChartNormal). n is chosen so
        that l is positive on every line going from the source to the receiver, hence on their positive combinations,
        the vertices of the cells: all of them land in the affine chart l = 1. The division is rational, without
        square root, and exact for the rational arithmetic.
      - a hyperplane (the line of an edge) can have l = 0, an edge parallel to the chart plane: it is scaled by the
        power of two that brings its largest coordinate to [0.5, 1), which is exact in floating point.*/
    template<class P, class S>
    class MathPluckerNormalization
    {
    public:
        MathPluckerNormalization()
            : mType(NO_NORMALIZATION), mChartNormal(MathVector3_<S>::Zero())
        {
        }

        MathPluckerNormalization(NormalizationType aType, const MathVector3_<S>& aChartNormal)
            : mType(aType), mChartNormal(aChartNormal)
        {
        }

        NormalizationType getType() const { return mType; }
        const MathVector3_<S>& getChartNormal() const { return mChartNormal; }

        /** @brief Whether the points are on the unit hypersphere, which the debug checks of the polyhedron verify*/
        bool isUnitNorm() const { return mType == HYPERSPHERE_NORMALIZATION; }

        /** @brief The value of the chart form at a point: its direction part along the chart normal*/
        S getChartValue(const P& aPoint) const
        {
            return aPoint.getDirection().dot(mChartNormal);
        }

        /** @brief Normalize a point of the cells: a vertex, or a line used as a point.

        In the chart, the point is divided by the absolute value of the chart form, which keeps the orientation of a
        line that would come reversed. A point whose form vanishes (a line parallel to the chart plane, which no line
        from the source to the receiver is) gets the power of two scale of the hyperplanes instead.*/
        P normalizePoint(const P& aPoint) const
        {
            if (mType == HYPERSPHERE_NORMALIZATION)
            {
                return aPoint.getNormalized();
            }
            return aPoint * getPointScale(aPoint);
        }

        /** @brief The positive factor normalizePoint multiplies a point by, for NONE and CHART (1 for NONE).*/
        S getPointScale(const P& aPoint) const
        {
            switch (mType)
            {
            case HYPERSPHERE_NORMALIZATION:
                return S(1) / MathArithmetic<S>::getSqrt(aPoint.getDirection().getSquaredNorm() + aPoint.getLocation().getSquaredNorm());
            case CHART_NORMALIZATION:
            {
                const S magnitude = MathArithmetic<S>::getAbs(getChartValue(aPoint));
                if (!(magnitude > MathArithmetic<S>::Tolerance()))
                {
                    return powerOfTwoScale(aPoint);
                }
                return S(1) / magnitude;
            }
            default:
                return S(1);
            }
        }

        /** @brief Normalize a symbolic point (see MathSoS): every term is multiplied by the factor of its leading term,
        which is one positive factor for the whole polynomial, the same projective point.*/
        void normalizeSymbolic(MathSymbolicPoint<P>& aPoint) const
        {
            if (mType == NO_NORMALIZATION)
            {
                return;
            }
            aPoint.scale(getPointScale(aPoint.leading()));
        }

        /** @brief Normalize a hyperplane, the line of an edge.*/
        P normalizeHyperplane(const P& aHyperplane) const
        {
            switch (mType)
            {
            case HYPERSPHERE_NORMALIZATION:
                return aHyperplane.getNormalized();
            case CHART_NORMALIZATION:
                return scaleByPowerOfTwo(aHyperplane);
            default:
                return aHyperplane;
            }
        }

        /** @brief The power of two scale bringing the largest absolute coordinate of a point to [0.5, 1).

        Multiplying by a power of two is exact in floating point: the point keeps its exact coordinates, only the
        exponent changes.*/
        static P scaleByPowerOfTwo(const P& aPoint)
        {
            return aPoint * powerOfTwoScale(aPoint);
        }

        /** @brief The power of two bringing the largest absolute coordinate of a point to [0.5, 1) (1 for a null point)*/
        static S powerOfTwoScale(const P& aPoint)
        {
            const MathVector3_<S>& d = aPoint.getDirection();
            const MathVector3_<S>& m = aPoint.getLocation();
            const double c[6] = { MathArithmetic<S>::to_double(d.x), MathArithmetic<S>::to_double(d.y), MathArithmetic<S>::to_double(d.z),
                                  MathArithmetic<S>::to_double(m.x), MathArithmetic<S>::to_double(m.y), MathArithmetic<S>::to_double(m.z) };
            double largest = 0.0;
            for (double v : c)
            {
                largest = std::fmax(largest, std::fabs(v));
            }
            if (!(largest > 0.0) || !std::isfinite(largest))
            {
                return S(1);
            }
            int exponent = 0;
            std::frexp(largest, &exponent);
            return S(std::ldexp(1.0, -exponent));
        }

        /** @brief The chart normal for a source and a receiver: a unit vector n such that (b - a) . n > 0 for every a of
        the source and b of the receiver.

        The query polygons are clipped so that each one lies on the positive side of the plane of the other. For a line
        from a point a of the source to a point b of the receiver, (b - a) . n0, with n0 the normal of the source plane,
        is the distance of b to that plane, positive; (b - a) . (-n1), with n1 the normal of the receiver plane, is the
        distance of a to it. The candidate kept is the one with the larger smallest distance, the better conditioned
        chart. aMargin returns that smallest distance: when it is not clearly positive, no chart is well conditioned
        and the caller falls back to the hypersphere.

        @param aSource: the source polygon, its plane oriented toward the receiver
        @param aReceiver: the receiver polygon, its plane oriented toward the source*/
        template<class POLYGON>
        static MathVector3_<S> computeChartNormal(const POLYGON& aSource, const POLYGON& aReceiver, double& aMargin)
        {
            const auto sourcePlane = aSource.getPlane();
            const auto receiverPlane = aReceiver.getPlane();

            double sourceMargin = std::numeric_limits<double>::max();
            for (int i = 0; i < aReceiver.getVertexCount(); i++)
            {
                sourceMargin = std::fmin(sourceMargin, sourcePlane.dot(aReceiver.getVertex(i)));
            }
            double receiverMargin = std::numeric_limits<double>::max();
            for (int i = 0; i < aSource.getVertexCount(); i++)
            {
                receiverMargin = std::fmin(receiverMargin, receiverPlane.dot(aSource.getVertex(i)));
            }

            // the plane normals are unit vectors: the dot products are distances
            if (sourceMargin >= receiverMargin)
            {
                aMargin = sourceMargin;
                return convert<MathVector3_<S>>(sourcePlane.mNormal);
            }
            aMargin = receiverMargin;
            return convert<MathVector3_<S>>(receiverPlane.mNormal * -1.0);
        }

    private:
        NormalizationType mType;
        MathVector3_<S> mChartNormal;
    };
}
