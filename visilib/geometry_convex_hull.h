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

#include <vector>
#include "math_vector_3.h"
#include "math_plane_3.h"
#include "visilib_core.h"
namespace visilib
{
    /** @brief Store a convex hull  under the form of a list of 3D plane equations*/

    class GeometryConvexHull
    {
    public:

        void addFace(const MathPlane3d& aFace)
        {
            mFaces.push_back(aFace);
        }
        const std::vector<MathPlane3d>& getFaces() const
        {
            return mFaces;
        }
    private:
        std::vector<MathPlane3d> mFaces;
    };

    /** @brief Compute the convex hull of a set of point in 3D space*/
    class GeometryConvexHullBuilder
    {
    public:
        static GeometryConvexHull* build(const std::vector<MathVector3d>& inputVertices0, const std::vector<MathVector3d>& inputVertices1);
    private:

     
        static bool compute(const std::vector<MathVector3d>& inputVertices0, const std::vector<MathVector3d>& inputVertices1, std::vector<MathPlane3d>& outputFaces);
    };

    inline GeometryConvexHull* GeometryConvexHullBuilder::build(const std::vector<MathVector3d>& inputVertices0, const std::vector<MathVector3d>& inputVertices1)
    {
        std::vector<MathPlane3d> myFaces;
        GeometryConvexHull* convexHull = nullptr;

        if (compute(inputVertices0, inputVertices1, myFaces)&& 
            compute(inputVertices1, inputVertices0, myFaces))
        {
            convexHull = new GeometryConvexHull();
            for (size_t i = 0; i < myFaces.size(); i++)
            {
                convexHull->addFace(myFaces[i]);
            }
        }

        return convexHull;

    }

    inline bool GeometryConvexHullBuilder::compute(const std::vector<MathVector3d>& inputVertices0, const std::vector<MathVector3d>& inputVertices1, std::vector<MathPlane3d>& outputFaces)
    {
        const size_t n0 = inputVertices0.size();
        const size_t n1 = inputVertices1.size();

        for (size_t i = 0; i < n0; i++)
        {
            if (n1>1)
            {
                for (size_t j = 0; j < n1; j++)
                {
                    {
                        outputFaces.push_back(MathGeometry::computePlane(inputVertices0[i], inputVertices1[(j+1)%n1], inputVertices1[j]));
                        const MathPlane3d & myFace = outputFaces[outputFaces.size() - 1];
                        for (size_t k = 0; k < n0; k++)
                        {
                            if (k != i)
                            {
                                double d = myFace.dot(inputVertices0[k]);
                                if (d < 0)
                                {
                                    outputFaces.pop_back();
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        return true;
    }
}
