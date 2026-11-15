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

#include "visilib_core.h"
#include "math_arithmetic.h"
#include "math_vector_3.h"
#include "math_plucker_traits.h"

namespace visilib
{

    /** @brief Represents a line parallel to a fixed direction, in Plucker space.

    The class models a set of parallel lines at large, the rays of a directional
    light: an instance is one of these rays, whose direction is given at
    construction and is never replaced afterwards, while the moment is free to
    move the ray along the screen orthogonal to it. The stored pair (d, b) is
    held to the equation d . b = 0, which is the Plucker quadric restricted to
    the family, and it holds it by construction: the projection of the moment
    orthogonal to the direction, applied by each operation that could break it,
    is the quadric correction of the generic space reduced to what the family
    needs. Every value of the class is therefore a real line, and nothing can be
    imaginary here.

    A silhouette edge is carried by the same class, the edge being a line with a
    direction of its own and a moment that is its true one: the Klein form taken
    between the edge and a ray of a family is the incidence of the two, the very
    linear half space test the ray is subject to, and the polyhedra and their
    predicates need no other operation on the class.

    @param S: the coordinates representation (float, double, exact)
    */
    template<class S>
    class MathPluckerParallel6
    {
    public:

        MathPluckerParallel6()
        {
        }

        /**@brief Construct the ray parallel to a given direction and going through a 3D point

        The moment p x d is orthogonal to the direction by construction, whatever
        the point: the entry into the family is exact for the rational arithmetics,
        and carries no square root.*/
        MathPluckerParallel6(const MathVector3_<S>& aDirection, const MathVector3_<S>& aPoint)
            :mDirection(aDirection),
            mLocation(MathVector3_<S>::cross(aPoint, aDirection))
        {
            V_ASSERT(MathArithmetic<S>::isFinite(getDirection()));
            V_ASSERT(MathArithmetic<S>::isFinite(getLocation()));
        }

        /**@brief Construct a ray of the family from another one*/
        MathPluckerParallel6(const MathPluckerParallel6<S>& v)
        {
            mDirection = v.mDirection;
            mLocation = v.mLocation;

            V_ASSERT(MathArithmetic<S>::isFinite(getDirection()));
            V_ASSERT(MathArithmetic<S>::isFinite(getLocation()));
        }

        /** @brief Carry a 3D edge, given by its two endpoints, as a plane of the family arithmetic

        The result holds the direction and the moment of the edge itself, the line
        going through the first endpoint with that direction being the edge: taken
        against a ray with the Klein form, it gives the incidence of the ray with
        the edge, which is the half space a silhouette imposes on the rays.*/
        static MathPluckerParallel6<S> throughEdge(const MathVector3_<S>& aBegin, const MathVector3_<S>& anEnd)
        {
            return MathPluckerParallel6<S>(anEnd - aBegin, aBegin);
        }

        /**@brief Add another line of the family in Plucker space

        The directions add up parallel to themselves and the sum of two moments
        orthogonal to the direction is orthogonal: the projection that follows is
        what keeps the result a line when an operand came from another family, and
        does nothing otherwise.*/
        MathPluckerParallel6<S>& operator -= (const MathPluckerParallel6<S>& v)
        {
            mDirection -= v.mDirection;
            mLocation -= v.mLocation;
            projectMoment();

            return *this;
        }

        /**@brief Add another line of the family in Plucker space*/
        MathPluckerParallel6<S>& operator += (const MathPluckerParallel6<S>& v)
        {
            mDirection += v.mDirection;
            mLocation += v.mLocation;
            projectMoment();

            return *this;
        }

        /**@brief Multiply all the coordinates by a provided number in Plucker space

        The line is unchanged, both parts scaling together: it is how a value of
        the family is renormalized without moving the ray.*/
        MathPluckerParallel6<S>& operator *= (S aNumber)
        {
            mLocation *= aNumber;
            mDirection *= aNumber;

            return *this;
        }

        /**@brief Return the addition of two lines of the family in Plucker space*/
        MathPluckerParallel6<S>  operator +(const MathPluckerParallel6<S>& v) const
        {
            MathPluckerParallel6<S>  result(*this);
            result += v;

            return result;
        }

        /**@brief Return the subtraction of two lines of the family in Plucker space*/
        MathPluckerParallel6<S>  operator -(const MathPluckerParallel6<S>& v) const
        {
            MathPluckerParallel6<S>  result(*this);
            result -= v;

            return result;
        }

        /**@brief Return a line with its coordinates divided by the given number in Plucker space*/
        MathPluckerParallel6<S>  operator / (S aNumber) const
        {
            MathPluckerParallel6<S>  result(*this);
            result.mLocation /= aNumber;
            result.mDirection /= aNumber;

            return result;
        }

        /**@brief Return a line with its coordinates multiplied by the given number in Plucker space*/
        MathPluckerParallel6<S>  operator * (S aNumber) const
        {
            MathPluckerParallel6<S>  result(*this);
            result.mLocation *= aNumber;
            result.mDirection *= aNumber;

            return result;
        }

        /**@brief Return the directional part of the Plucker point

        It is the direction of the family, parallel to what it was given at
        construction, up to the scale the arithmetic left on it.*/
        const MathVector3_<S>& getDirection()const
        {
            return mDirection;
        }


        /**@brief Return the location part of the Plucker point

        It is the moment of the line, orthogonal to the direction by the invariant
        of the class, and a coordinate of the ray on the screen of the family.*/
        const MathVector3_<S>& getLocation()const
        {
            return mLocation;
        }

        /** @brief Compute the Plucker dot product of the line with another line of the space

        The sign of the result provides the relative orientation of the two lines
        (0: intersect, <0: skew ccw, >0: skew cw). Against an edge carried by
        throughEdge, it is the incidence of the ray with the edge, which is what
        decides on which side of a silhouette plane a line of the cell lies.
        */
        S dot(const MathPluckerParallel6<S>& p2) const
        {
            return mDirection.dot(p2.mLocation) + mLocation.dot(p2.mDirection);
        }

        /** @brief Compute the ordinary scalar product of the two lines in the six dimensional Plucker space

        Not the Klein form that dot gives, which vanishes for a line with itself:
        taken with itself it is the squared six dimensional norm, the scale that
        renormalization divides by.*/
        S scalarDot(const MathPluckerParallel6<S>& p2) const
        {
            return mDirection.dot(p2.mDirection) + mLocation.dot(p2.mLocation);
        }

        /** @brief Return the value of the Plucker quadric at this line, which is nought

        The family lies in the quadric: d . b = 0 is the equation of the quadric
        restricted to lines parallel to d, and the class holds its moments to it.
        The method exists so that a generic predicate of the six dimensional space
        can be asked here, and answers that there is nothing to correct.*/
        S getQuadric() const
        {
            return mDirection.dot(mLocation);
        }

        /** @brief Compute the gradient of the equation of the Plucker quadric at this line

        The quadric is the zero set of Q(a, b) = a . b, whose gradient at (a, b)
        is (b, a). The two parts swap, and a moment orthogonal to its direction is
        a direction orthogonal to its moment: the gradient of a line of a family
        is a line of the family orthogonal to it, and satisfies the invariant.
        */
        MathPluckerParallel6<S> quadricGradient() const
        {
            return MathPluckerParallel6<S>(mLocation, mDirection, RawPair());
        }

        /** @brief Return a normalized line by projection on the unit hypersphere

        Both parts scale together, the line is unchanged and so is its screen
        point: this is a renormalization of the coordinates, not a move.*/
        MathPluckerParallel6<S> getNormalized() const
        {
            V_ASSERT(MathArithmetic<S>::isFinite(getDirection()));
            V_ASSERT(MathArithmetic<S>::isFinite(getLocation()));

            MathPluckerParallel6<S>  result(*this);

            S myDiv = MathArithmetic<S>::getSqrt(mDirection.getSquaredNorm() + mLocation.getSquaredNorm());
            V_ASSERT(MathArithmetic<S>::getAbs(myDiv) > MathArithmetic<S>::Tolerance());

            result.mDirection /= myDiv;
            result.mLocation /= myDiv;

            V_ASSERT(MathArithmetic<S>::isFinite(result.getDirection()));
            V_ASSERT(MathArithmetic<S>::isFinite(result.getLocation()));

            return result;
        }

        /** @brief Return the point where the line meets the screen of the family

        The screen is the plane through the origin orthogonal to the direction,
        and the moment is its coordinate: o = (d x b) / (d . d) recovers the 3D
        point of the ray, p x d giving back the moment from any point of it.*/
        MathVector3_<S> getScreenPoint() const
        {
            S d2 = mDirection.getSquaredNorm();
            V_ASSERT(d2 > MathArithmetic<S>::Tolerance());

            MathVector3_<S> result = MathVector3_<S>::cross(mDirection, mLocation);
            result /= d2;

            return result;
        }

        /** @brief Return the distance between two lines of the same family

        It is the offset of their screen points, the rays never meeting; between
        lines of different families, the generic formula of the six dimensional
        space, the skew distance the Klein form gives.*/
        S getDistance(const MathPluckerParallel6<S>& other) const
        {
            S myDiv = MathVector3_<S>::cross(getDirection(), other.getDirection()).getNorm();

            if (MathArithmetic<S>::getAbs(myDiv) > MathArithmetic<S>::Tolerance())
            {
                return MathArithmetic<S>::getAbs(dot(other)) / myDiv;
            }
            else
            {
                return (getScreenPoint() - other.getScreenPoint()).getNorm();
            }
        }

        /**@brief The line through the origin of the family, which is also the null object

        The ray parallel to a nought direction is none, but the arithmetic needs a
        zero to start a sum from, and adding to it yields the other operand.*/
        static MathPluckerParallel6<S>  Zero()
        {
            return MathPluckerParallel6<S>();
        }

        bool isZero(S anEps)const
        {
            if (!getLocation().isZero(anEps))
                return false;
            if (!getDirection().isZero(anEps))
                return false;
            return true;
        }

        friend std::ostream& operator<<(std::ostream& os, const MathPluckerParallel6& v)
        {
            os << v.mDirection << v.mLocation;
            return os;
        }

        MathPluckerParallel6<S>& operator = (const MathPluckerParallel6<S>& v)
        {
            mDirection = v.mDirection;
            mLocation = v.mLocation;

            return *this;
        }

    private:

        /* The two parts of a line, adopted as they are rather than rebuilt from a
        point: the gradient of the quadric swaps them, and a swapped pair is a
        direction and a moment no ray is walked through.*/
        struct RawPair {};

        MathPluckerParallel6(const MathVector3_<S>& aDirection, const MathVector3_<S>& aMoment, RawPair)
            : mDirection(aDirection),
            mLocation(aMoment)
        {
            V_ASSERT(MathArithmetic<S>::isFinite(getDirection()));
            V_ASSERT(MathArithmetic<S>::isFinite(getLocation()));
        }

        /* Bring the moment back orthogonal to the direction: what a six dimensional
        sum leaves stretched along the direction is not a line any more, and this is
        the whole of the quadric correction the family knows, a projection and no
        square root.*/
        void projectMoment()
        {
            const S d2 = mDirection.getSquaredNorm();
            if (d2 <= MathArithmetic<S>::Tolerance())
                return;

            mLocation -= mDirection * (mDirection.dot(mLocation) / d2);
        }

        MathVector3_<S> mDirection;   /**< @brief Direction part of the Plucker point, fixed at construction*/
        MathVector3_<S> mLocation;    /**< @brief Moment of the line, the coordinate moving the ray on the screen*/

    };

    typedef MathPluckerParallel6<float> MathPluckerParallel6f;
    typedef MathPluckerParallel6<double> MathPluckerParallel6d;

    /* The family is the case the traits header describes: every one of its lines is
    real, and its cells are polygons on the screen, vertices at two lines, edges
    joining vertices that share one. */

    template<class S>
    struct MathPluckerCellTraits<MathPluckerParallel6<S>>
    {
        static constexpr size_t cCommonFacetsForEdge = 1;
        static constexpr bool cAllLinesReal = true;
    };

    /* A 3D edge is not a ray of the family, its direction is its own: the edge comes
    in as the carrier throughEdge describes it, the pair the Klein form reads against
    a ray, which is the half space the edge imposes on the screen. */

    template<class S>
    struct PluckerEdgeCarrier<MathPluckerParallel6<S>>
    {
        static MathPluckerParallel6<S> fromPoints(const MathVector3d& a, const MathVector3d& b)
        {
            return MathPluckerParallel6<S>::throughEdge(convert<MathVector3_<S>>(a), convert<MathVector3_<S>>(b));
        }
    };
}
