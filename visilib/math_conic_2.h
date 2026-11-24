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

#include "math_arithmetic.h"
#include "math_frame_3.h"
#include "math_plucker_6.h"
#include "math_vector_2.h"
#include "math_vector_3.h"

namespace visilib
{
    template<class S>
    class MathConic2_;

    typedef MathConic2_<float> MathConic2f;
    typedef MathConic2_<double> MathConic2d;

    /** @brief A conic of the plane, the zero set of the second order polynomial

        A u^2 + B u v + C v^2 + D u + E v + F = 0

    The coefficients are defined up to a common factor.
        @param S the coordinates representation (float, double, exact)*/
    template<class S>
    class MathConic2_
    {
    public:
        enum Coefficient { A = 0, B, C, D, E, F, COUNT };

        MathConic2_()
        {
            for (int i = 0; i < COUNT; i++)
            {
                mCoefficients[i] = S(0);
            }
        }

        MathConic2_(S a, S b, S c, S d, S e, S f)
        {
            mCoefficients[A] = a; mCoefficients[B] = b; mCoefficients[C] = c;
            mCoefficients[D] = d; mCoefficients[E] = e; mCoefficients[F] = f;
        }

        /** @brief The value of the polynomial at aPoint: zero on the conic, its sign tells the two sides apart.*/
        S evaluate(const MathVector2_<S>& aPoint) const
        {
            const S u = aPoint.x;
            const S v = aPoint.y;
            return mCoefficients[A] * u * u + mCoefficients[B] * u * v + mCoefficients[C] * v * v
                + mCoefficients[D] * u + mCoefficients[E] * v + mCoefficients[F];
        }

        /** @brief The gradient of the polynomial at aPoint, the normal direction of the conic there.*/
        MathVector2_<S> getGradient(const MathVector2_<S>& aPoint) const
        {
            const S u = aPoint.x;
            const S v = aPoint.y;
            return MathVector2_<S>(2 * mCoefficients[A] * u + mCoefficients[B] * v + mCoefficients[D],
                                   mCoefficients[B] * u + 2 * mCoefficients[C] * v + mCoefficients[E]);
        }

        /** @brief The distance of aPoint to the conic, to first order: |f| / |grad f|. Infinite where the gradient vanishes.*/
        S getAlgebraicDistance(const MathVector2_<S>& aPoint) const
        {
            const S g = getGradient(aPoint).getNorm();
            const S f = MathArithmetic<S>::getAbs(evaluate(aPoint));
            if (g <= S(0))
            {
                return f <= S(0) ? S(0) : S(std::numeric_limits<double>::infinity());
            }
            return f / g;
        }

        /** @brief Divide the coefficients by the greatest of their absolute values. Return that value (zero leaves the conic as is).*/
        S normalize()
        {
            S scale = S(0);
            for (int i = 0; i < COUNT; i++)
            {
                const S c = MathArithmetic<S>::getAbs(mCoefficients[i]);
                if (c > scale)
                {
                    scale = c;
                }
            }
            if (scale > S(0))
            {
                for (int i = 0; i < COUNT; i++)
                {
                    mCoefficients[i] /= scale;
                }
            }
            return scale;
        }

        /** @brief Whether the second order part vanishes: the conic is then a line (or empty), D u + E v + F = 0.*/
        bool isLinear(const S& tolerance) const
        {
            return MathArithmetic<S>::getAbs(mCoefficients[A]) <= tolerance
                && MathArithmetic<S>::getAbs(mCoefficients[B]) <= tolerance
                && MathArithmetic<S>::getAbs(mCoefficients[C]) <= tolerance;
        }

        const S& operator[](int i) const { return mCoefficients[i]; }
        S& operator[](int i) { return mCoefficients[i]; }

        /** @brief The trace, on the plane of aFrame, of the surface swept by the transversals of three lines.

        A point x lies on a line meeting the three lines (d_i, m_i) when some direction t is orthogonal to the three
        vectors (x - a_i) x d_i = x x d_i - m_i, a_i any point of line i: when their determinant vanishes,

            F(x) = det [ x x d_i - m_i ] = 0

        The terms of third order cancel ((x x d_2) x (x x d_3) is along x, orthogonal to x x d_1), so F is of second
        order: for three skew lines it is the regulus they span, for lines in a degenerate position a pair of planes.
        On the plane x = O + u e_u + v e_v, each row is affine in (u, v): row_i = r_i0 + u r_iu + v r_iv, and the
        determinant expands, by multilinearity, over the 27 choices of one term per row. The terms of third order
        are dropped; the conic is returned normalized (see normalize).*/
        static MathConic2_ fromTransversalsOfThreeLines(const MathPlucker6<S> lines[3], const MathFrame3_<S>& aFrame)
        {
            MathVector3_<S> rows[3][3]; // [line][0: constant, 1: u, 2: v]
            for (int i = 0; i < 3; i++)
            {
                const MathVector3_<S>& d = lines[i].getDirection();
                rows[i][0] = MathVector3_<S>::cross(aFrame.getOrigin(), d) - lines[i].getLocation();
                rows[i][1] = MathVector3_<S>::cross(aFrame.getAxisU(), d);
                rows[i][2] = MathVector3_<S>::cross(aFrame.getAxisV(), d);
            }

            // monomial[degree of u][degree of v]
            S monomial[4][4];
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    monomial[i][j] = S(0);
                }
            }

            for (int c0 = 0; c0 < 3; c0++)
            {
                for (int c1 = 0; c1 < 3; c1++)
                {
                    for (int c2 = 0; c2 < 3; c2++)
                    {
                        const int choice[3] = { c0, c1, c2 };
                        int du = 0;
                        int dv = 0;
                        for (int k = 0; k < 3; k++)
                        {
                            du += choice[k] == 1 ? 1 : 0;
                            dv += choice[k] == 2 ? 1 : 0;
                        }
                        if (du + dv > 2)
                        {
                            continue;
                        }
                        const S det = rows[0][c0].dot(MathVector3_<S>::cross(rows[1][c1], rows[2][c2]));
                        monomial[du][dv] += det;
                    }
                }
            }

            MathConic2_ conic(monomial[2][0], monomial[1][1], monomial[0][2], monomial[1][0], monomial[0][1], monomial[0][0]);
            conic.normalize();
            return conic;
        }

    private:
        S mCoefficients[COUNT];
    };
}
