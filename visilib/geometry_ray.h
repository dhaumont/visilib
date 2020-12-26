/*
Visilib, an open source library for exact visibility computation.
Copyright(C) 2020 by Denis Haumont

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

#include "math_vector_3.h"
#include "visilib.h"

namespace visilib
{
    /** @brief Store a ray in 3D space (ie a vertex and a direction). It also stores the coefficients required for watertight ray-triangle intersection */

    class GeometryRay
    {
    public:

        GeometryRay(const Ray& ray);

        GeometryRay(const MathVector3f& aStart, const MathVector3f& aDirection);

        GeometryRay(const GeometryRay& other);

        void initialize(const MathVector3f& aStart, const MathVector3f& aDirection);

        /** @brief Get the starting point of the ray*/
        const MathVector3f& getStart() const
        {
            return mStart;
        }

        /** @brief Get the direction of the ray*/
        const MathVector3f& getDirection() const
        {
            return mDirection;
        }

        /** @brief Get the coefficients for watertight ray-triangle intersection (double precision)*/
        void get(double& aSx, double& aSy, double& aSz, int& akx, int& aky, int& akz) const
        {
            aSx = Sx;
            aSy = Sy;
            aSz = Sz;
            akx = kx;
            aky = ky;
            akz = kz;
        }

        /** @brief Get the coefficients for watertight ray-triangle intersection (single precision)*/
        void get(float& aSx, float& aSy, float& aSz, int& akx, int& aky, int& akz) const
        {
            aSx = (float)Sx;
            aSy = (float)Sy;
            aSz = (float)Sz;
            akx = kx;
            aky = ky;
            akz = kz;
        }
    private:

        /** @brief The starting point of the ray*/
        MathVector3f mStart;

        /** @brief The direction of the ray */
        MathVector3f mDirection;

        /** @brief The coefficients for watertight ray-triangle intersection */
        double Sx;
        double Sy;
        double Sz;
        int kx;
        int ky;
        int kz;
    };

    inline GeometryRay::GeometryRay(const Ray& ray)
        : Sx(0.0), Sy(0.0), Sz(0.0), kx(0), ky(0), kz(0)
    {
        initialize(MathVector3f(ray.org), MathVector3f(ray.dir));
    }

    inline GeometryRay::GeometryRay(const MathVector3f& aStart, const MathVector3f& aDirection)
        : Sx(0.0), Sy(0.0), Sz(0.0), kx(0), ky(0), kz(0)
    {
        initialize(aStart, aDirection);
    }

    inline GeometryRay::GeometryRay(const GeometryRay& other)
        : Sx(0.0), Sy(0.0), Sz(0.0), kx(0), ky(0), kz(0)
    {
        initialize(other.mStart, other.mDirection);
    }

    // Note: the initialize method is defined in "MathGeometry.h"
}