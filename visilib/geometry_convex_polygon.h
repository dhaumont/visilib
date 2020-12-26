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

#include <vector>
#include "math_vector_2.h"
#include "math_vector_3.h"
#include "math_plane_3.h"

namespace visilib
{
    /** @brief Store an convex polygon in 3D space*/

    class GeometryConvexPolygon
    {
    public:
        GeometryConvexPolygon(const MathVector3f& v0, const MathVector3f& v1, const MathVector3f& v2);
        GeometryConvexPolygon(const std::vector<MathVector3d>& aList);
        GeometryConvexPolygon(const std::vector<MathVector3f>& aList);
        GeometryConvexPolygon(const GeometryConvexPolygon& aPoly);
        GeometryConvexPolygon(MathVector3d* anArray, size_t aNumber);
        GeometryConvexPolygon(const float* anArray, size_t numVertices);

        MathVector2i getEdge(size_t i) const
        {
            return MathVector2i(i, (i + 1) % getVertexCount());
        }

        MathVector2i getEdgesOfVertex(size_t aVertexNumber) const
        {
            if (aVertexNumber != 0)
            {
                return MathVector2i(aVertexNumber, aVertexNumber - 1);
            }
            else
            {
                return MathVector2i(getVertexCount() - 1, 0);
            }
        }

        const std::vector<MathVector3d>& getVertices() const
        {
            return mVertices;
        }

        const MathVector3d& getVertex(size_t i)const;

        const MathPlane3d& getPlane() const
        {
            return mPlane;
        }

        void setPlane(const MathPlane3d & aPlane)
        {
            mPlane = aPlane;
        }

        size_t getVertexCount() const
        {
            return mVertices.size();
        }

        void clear()
        {
            mVertices.clear();
            mPlane = MathPlane3d::Zero();
        }

        bool isValid()const;

    private:
        void computePlane();

        std::vector<MathVector3d> mVertices;
        MathPlane3d mPlane;
    };

    inline GeometryConvexPolygon::GeometryConvexPolygon(const MathVector3f & v0, const MathVector3f & v1, const MathVector3f & v2)
    {
        mVertices.push_back(convert<MathVector3d>(v0));
        mVertices.push_back(convert<MathVector3d>(v1));
        mVertices.push_back(convert<MathVector3d>(v2));
        computePlane();
    }

    inline GeometryConvexPolygon::GeometryConvexPolygon(const GeometryConvexPolygon & other)
    {
        mVertices = other.mVertices;
        mPlane = other.mPlane;
    }

    inline GeometryConvexPolygon::GeometryConvexPolygon(const std::vector<MathVector3d> & aList)
    {
        mVertices = aList;
        computePlane();
    }

    inline GeometryConvexPolygon::GeometryConvexPolygon(const std::vector<MathVector3f> & aList)
    {
        for (size_t i = 0; i < aList.size(); i++)
        {
            mVertices.push_back(convert<MathVector3d>(aList[i]));
        }

        computePlane();
    }

    inline GeometryConvexPolygon::GeometryConvexPolygon(const float* anArray, size_t numVertices)
    {
        for (size_t i = 0; i < numVertices; i++)
        {
            mVertices.push_back(MathVector3d(anArray[i * 3 + 0], anArray[i * 3 + 1], anArray[i * 3 + 2]));
        }

        computePlane();
    }

    inline const MathVector3d& GeometryConvexPolygon::getVertex(size_t i)const
    {
        return mVertices[i];
    }

    inline bool GeometryConvexPolygon::isValid()const
    {
        return getVertexCount() > 0;
    }
}