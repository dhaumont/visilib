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
#include "math_plane_3.h"

namespace visilib
{

    /** @brief Represents a point in 6D Plucker space.

    The 6D points in Plucker space are encoded as two 3D vectors.
        @param S: the coordinates representation (float, double, exact)
    */
    template<class S>
    class MathPlucker6
    {
    public:

        MathPlucker6()
        {
        }

        /**@brief Construct a plucker point from its 6 coefficients in Plucker space*/
        MathPlucker6(S x0, S x1, S x2, S x3, S x4, S x5)
            :mDirection(x0, x1, x2),
            mLocation(x3, x4, x5)
        {
            V_ASSERT(MathArithmetic<S>::isFinite(getDirection()));
            V_ASSERT(MathArithmetic<S>::isFinite(getLocation()));
        }

        /**@brief Construct a plucker point from another point*/
        MathPlucker6(const MathPlucker6<S>& v)
        {
            mDirection = v.mDirection;
            mLocation = v.mLocation;

            V_ASSERT(MathArithmetic<S>::isFinite(getDirection()));
            V_ASSERT(MathArithmetic<S>::isFinite(getLocation()));
        }

        /**@brief Construct a plucker point from a 3D line defined by two points*/
        MathPlucker6(const MathVector3d& aBegin, const MathVector3d& anEnd)
        {
            MathVector3_<S> myBegin = convert<MathVector3_<S> >(aBegin);
            MathVector3_<S> myEnd = convert<MathVector3_<S> >(anEnd);

            mDirection = myEnd - myBegin;
            mLocation = MathVector3_<S>::cross(myBegin, myEnd);

            V_ASSERT(MathArithmetic<S>::isFinite(getDirection()));
            V_ASSERT(MathArithmetic<S>::isFinite(getLocation()));
        }

        /**@brief Construct a plucker point from a 3D line defined by two points*/
        MathPlucker6(const MathVector3f& aBegin, const MathVector3f& anEnd)
        {
            MathVector3_<S> myBegin = convert<MathVector3_<S>>(aBegin);
            MathVector3_<S> myEnd = convert<MathVector3_<S>>(anEnd);

            mDirection = myEnd - myBegin;
            mLocation = MathVector3_<S>::cross(myBegin, myEnd);

            V_ASSERT(MathArithmetic<S>::isFinite(getDirection()));
            V_ASSERT(MathArithmetic<S>::isFinite(getLocation()));
        }

        /**@brief Subtract another point in Plucker space*/
        MathPlucker6<S>& operator -= (const MathPlucker6<S>& v)
        {
            mLocation -= v.mLocation;
            mDirection -= v.mDirection;

            return *this;
        }

        /**@brief Add another point in Plucker space*/
        MathPlucker6<S>& operator += (const MathPlucker6<S>& v)
        {
            mLocation += v.mLocation;
            mDirection += v.mDirection;

            return *this;
        }

        /**@brief Multiply all the coordinates by a provided number in Plucker space*/
        MathPlucker6<S>& operator *= (S aNumber)
        {
            mLocation *= aNumber;
            mDirection *= aNumber;

            return *this;
        }

        /**@brief Return the addition of two points in Plucker space*/
        MathPlucker6<S>  operator +(const MathPlucker6<S>& v) const
        {
            MathPlucker6<S>  result(*this);
            result += v;

            return result;
        }

        /**@brief Return the subtraction of two points in Plucker space*/
        MathPlucker6<S>  operator -(const MathPlucker6<S>& v) const
        {
            MathPlucker6<S>  result(*this);
            result -= v;

            return result;
        }

        /**@brief Return a point with its coordinates divided by the given number in Plucker space*/
        MathPlucker6<S>  operator / (S aNumber) const
        {
            MathPlucker6<S>  result(*this);
            result.mLocation /= aNumber;
            result.mDirection /= aNumber;

            return result;
        }

        /**@brief Rerturn a point with its coordinates multiplied by the given number in Plucker space*/
        MathPlucker6<S>  operator * (S aNumber) const
        {
            MathPlucker6<S>  result(*this);
            result.mLocation *= aNumber;
            result.mDirection *= aNumber;

            return result;
        }

        /**@brief Return the directional part of the Plucker point*/
        const MathVector3_<S>& getDirection()const
        {
            return mDirection;
        }


        /**@brief Return the location part of the Plucker point*/
        const MathVector3_<S>& getLocation()const
        {
            return mLocation;
        }

        /** @brief Compute the Plucker dot product of the point with another point in Plucker space

        The sign of the result provides the relative orientation of the two lines (0: intersect, <0: skew ccw, >0: skew cw)
        */
        S dot(const MathPlucker6<S>& p2) const
        {
            return mDirection.dot(p2.mLocation) + mLocation.dot(p2.mDirection);
        }

        /** @brief The line through aPoint with the direction aDirection (the location part is aPoint x aDirection).*/
        static MathPlucker6<S> fromPointAndDirection(const MathVector3_<S>& aPoint, const MathVector3_<S>& aDirection)
        {
            MathPlucker6<S> result;
            result.mDirection = aDirection;
            result.mLocation = MathVector3_<S>::cross(aPoint, aDirection);
            return result;
        }

        /** @brief A point of the line: the foot of the perpendicular from the origin, (d x m) / |d|^2. The line must have a direction.*/
        MathVector3_<S> getClosestPointToOrigin() const
        {
            MathVector3_<S> p = MathVector3_<S>::cross(mDirection, mLocation);
            p /= mDirection.getSquaredNorm();
            return p;
        }

        /** @brief Whether the two lines meet at a finite point.

        The Klein form (dot) vanishes for coplanar lines. Parallel lines are coplanar but meet at infinity: their
        direction cross product vanishes as well, and they are reported as not meeting. The Klein form is compared
        relatively, to the product of the direction and location norms, so that the test does not depend on the scale
        of the coordinates.*/
        bool meets(const MathPlucker6<S>& other, const S& tolerance) const
        {
            const MathVector3_<S> n = MathVector3_<S>::cross(mDirection, other.mDirection);
            const S scale = mDirection.getSquaredNorm() * other.mDirection.getSquaredNorm();
            if (n.getSquaredNorm() <= tolerance * scale)
            {
                return false;
            }
            const S klein = dot(other);
            const S kleinScale = mDirection.getNorm() * other.mLocation.getNorm() + mLocation.getNorm() * other.mDirection.getNorm() + n.getNorm();
            return MathArithmetic<S>::getAbs(klein) <= tolerance * kleinScale;
        }

        /** @brief The point where the two lines meet, or false when they are parallel.

        With p a point of this line, the meeting point p + lambda d satisfies (p + lambda d) x d' = m', so that
        lambda (d x d') = m' - p x d', solved in the least square sense along n = d x d'. The point is only meaningful
        when the lines meet (see meets): for skew lines it is the point of this line closest to the other one along n.*/
        bool getMeetingPoint(const MathPlucker6<S>& other, MathVector3_<S>& aPoint, const S& tolerance) const
        {
            const MathVector3_<S> n = MathVector3_<S>::cross(mDirection, other.mDirection);
            const S nn = n.getSquaredNorm();
            if (nn <= tolerance * mDirection.getSquaredNorm() * other.mDirection.getSquaredNorm())
            {
                return false;
            }
            const MathVector3_<S> p = getClosestPointToOrigin();
            const S lambda = (other.mLocation - MathVector3_<S>::cross(p, other.mDirection)).dot(n) / nn;
            aPoint = p + mDirection * lambda;
            return true;
        }

        /** @brief Compute the ordinary scalar product of the two points in the six dimensional Plucker space

        This is not the Klein form that dot gives, which vanishes exactly for the pairs of lines meeting: taken
        with itself it is the squared six dimensional norm, which is what the projections onto the quadric and
        onto the sphere are made of.*/
        S scalarDot(const MathPlucker6<S>& p2) const
        {
            return mDirection.dot(p2.mDirection) + mLocation.dot(p2.mLocation);
        }

        /** @brief Compute the gradient of the equation of the Plucker quadric at this point

        The quadric is the zero set of Q(a, b) = a . b, whose gradient at (a, b) is (b, a).
        */
        MathPlucker6<S> quadricGradient() const
        {
            return MathPlucker6<S>(mLocation.x, mLocation.y, mLocation.z, mDirection.x, mDirection.y, mDirection.z);
        }

        /** @brief Return a normalized Plucker point by projection on the unit hypersphere */
        MathPlucker6<S> getNormalized() const
        {
            V_ASSERT(MathArithmetic<S>::isFinite(getDirection()));
            V_ASSERT(MathArithmetic<S>::isFinite(getLocation()));

            MathPlucker6<S>  result(*this);

            S myDiv = MathArithmetic<S>::getSqrt(mDirection.getSquaredNorm() + mLocation.getSquaredNorm());
            V_ASSERT(MathArithmetic<S>::getAbs(myDiv) > MathArithmetic<S>::Tolerance());

            result.mDirection /= myDiv;
            result.mLocation /= myDiv;

            V_ASSERT(MathArithmetic<S>::isFinite(result.getDirection()));
            V_ASSERT(MathArithmetic<S>::isFinite(result.getLocation()));

            return result;
        }

        S getDistance(const MathPlucker6<S>& other) const
        {
            S myDiv = MathVector3_<S>::cross(getDirection(),other.getDirection()).getNorm();

            if (MathArithmetic<S>::getAbs(myDiv) > MathArithmetic<S>::Tolerance())
            {
                return MathArithmetic<S>::getAbs(dot(other)) / myDiv;
            }
            else
            {
                // parallel
                S l1l1 = 1.0 / getDirection().getSquaredNorm();
                S l2 = other.getDirection().getNorm();

                S s = l2 * sqrt(l1l1);
                MathVector3_<S> l = other.getLocation();
                l /= s;
                return MathVector3_<S>::cross(getDirection(),getLocation() - l).getNorm() * l1l1;
            }
        }


        /**@brief Origin point in Plucker space*/
        static MathPlucker6<S>  Zero()
        {
            return MathPlucker6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        }

        bool isZero(S anEps)const
        {
           if (!getLocation().isZero(anEps))
              return false;
           if (!getDirection().isZero(anEps))
              return false;
           return true;
        }

        friend std::ostream& operator<<(std::ostream& os, const MathPlucker6& v)
        {
            os << v.mDirection << v.mLocation;
            return os;
        }

    private:
        MathVector3_<S> mDirection;   /**< @brief Direction part of the Plucker point*/
        MathVector3_<S> mLocation;    /**< @brief Location part of the Plucker point*/

    };

    typedef MathPlucker6<float> MathPlucker6f;
    typedef MathPlucker6<double> MathPlucker6d;

    template<>
    inline MathPlucker6f convert(const MathPlucker6d& f)
    {
        const MathVector3d& a = f.getDirection();
        const MathVector3d& b = f.getLocation();

        return MathPlucker6f((float)a.x, (float)a.y, (float)a.z, (float)b.x, (float)b.y, (float)b.z);
    }

    template<>
    inline MathPlucker6d convert(const MathPlucker6f& f)
    {
        const MathVector3f& a = f.getDirection();
        const MathVector3f& b = f.getLocation();

        return MathPlucker6d((double)a.x, (double)a.y, (double)a.z, (double)b.x, (double)b.y, (double)b.z);
    }

    template<>
    inline MathPlucker6f convert(const MathPlucker6f& f)
    {
        return f;
    }

    template<>
    inline MathPlucker6d convert(const MathPlucker6d& f)
    {
        return f;
    }

#ifdef ENABLE_LEDA
    template<>
    inline MathPlucker6d convert(const MathPlucker6<MathLedaReal>& f)
    {
        const MathVector3_<MathLedaReal>& a = f.getDirection();
        const MathVector3_<MathLedaReal>& b = f.getLocation();

        return MathPlucker6d(
            MathArithmetic<MathLedaReal>::to_double(a.x),
            MathArithmetic<MathLedaReal>::to_double(a.y),
            MathArithmetic<MathLedaReal>::to_double(a.z),
            MathArithmetic<MathLedaReal>::to_double(b.x),
            MathArithmetic<MathLedaReal>::to_double(b.y),
            MathArithmetic<MathLedaReal>::to_double(b.z));
    }

    template<>
    inline MathPlucker6<MathLedaReal> convert(const MathPlucker6d& f)
    {
        const MathVector3d& a = f.getDirection();
        const MathVector3d& b = f.getLocation();

        return MathPlucker6<MathLedaReal>(a.x, a.y, a.z, b.x, b.y, b.z);
    }
#endif

#if ENABLE_GMP
    template<>
    inline MathPlucker6d convert(const MathPlucker6<MathGmpFloat>& f)
    {
        const MathVector3_<MathGmpFloat>& a = f.getDirection();
        const MathVector3_<MathGmpFloat>& b = f.getLocation();

        return MathPlucker6d(
            MathArithmetic<MathGmpFloat>::to_double(a.x),
            MathArithmetic<MathGmpFloat>::to_double(a.y),
            MathArithmetic<MathGmpFloat>::to_double(a.z),
            MathArithmetic<MathGmpFloat>::to_double(b.x),
            MathArithmetic<MathGmpFloat>::to_double(b.y),
            MathArithmetic<MathGmpFloat>::to_double(b.z));
    }

    template<>
    inline MathPlucker6<MathGmpFloat> convert(const MathPlucker6d& f)
    {
        const MathVector3d& a = f.getDirection();
        const MathVector3d& b = f.getLocation();

        return MathPlucker6<MathGmpFloat>(a.x, a.y, a.z, b.x, b.y, b.z);
    }

    template<>
    inline MathPlucker6d convert(const MathPlucker6<MathGmpRational>& f)
    {
        const MathVector3_<MathGmpRational>& a = f.getDirection();
        const MathVector3_<MathGmpRational>& b = f.getLocation();

        return MathPlucker6d(
            MathArithmetic<MathGmpRational>::to_double(a.x),
            MathArithmetic<MathGmpRational>::to_double(a.y),
            MathArithmetic<MathGmpRational>::to_double(a.z),
            MathArithmetic<MathGmpRational>::to_double(b.x),
            MathArithmetic<MathGmpRational>::to_double(b.y),
            MathArithmetic<MathGmpRational>::to_double(b.z));
    }

    template<>
    inline MathPlucker6<MathGmpRational> convert(const MathPlucker6d& f)
    {
        const MathVector3d& a = f.getDirection();
        const MathVector3d& b = f.getLocation();

        return MathPlucker6<MathGmpRational>(a.x, a.y, a.z, b.x, b.y, b.z);
    }
#endif
#if defined(ENABLE_GMP) && defined(ENABLE_MPFI)
    template<>
    inline MathPlucker6d convert(const MathPlucker6<MathLazyRational>& f)
    {
        const MathVector3_<MathLazyRational>& a = f.getDirection();
        const MathVector3_<MathLazyRational>& b = f.getLocation();

        return MathPlucker6d(a.x.to_double(), a.y.to_double(), a.z.to_double(), b.x.to_double(), b.y.to_double(), b.z.to_double());
    }

    template<>
    inline MathPlucker6<MathLazyRational> convert(const MathPlucker6d& f)
    {
        const MathVector3d& a = f.getDirection();
        const MathVector3d& b = f.getLocation();

        return MathPlucker6<MathLazyRational>(a.x, a.y, a.z, b.x, b.y, b.z);
    }
#endif

#if ENABLE_MPFR
    template<>
    inline MathPlucker6d convert(const MathPlucker6<MathMpfr>& f)
    {
        const MathVector3_<MathMpfr>& a = f.getDirection();
        const MathVector3_<MathMpfr>& b = f.getLocation();

        return MathPlucker6d(
            MathArithmetic<MathMpfr>::to_double(a.x),
            MathArithmetic<MathMpfr>::to_double(a.y),
            MathArithmetic<MathMpfr>::to_double(a.z),
            MathArithmetic<MathMpfr>::to_double(b.x),
            MathArithmetic<MathMpfr>::to_double(b.y),
            MathArithmetic<MathMpfr>::to_double(b.z));
    }

    template<>
    inline MathPlucker6<MathMpfr> convert(const MathPlucker6d& f)
    {
        const MathVector3d& a = f.getDirection();
        const MathVector3d& b = f.getLocation();

        return MathPlucker6<MathMpfr>(a.x, a.y, a.z, b.x, b.y, b.z);
    }
#endif
}
