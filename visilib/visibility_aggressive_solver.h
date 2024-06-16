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

#include "visibility_solver.h"
#include "visibility_exact_query.h"

namespace visilib
{
    /** @brief TODO */

    template<class P, class S>
    class VisibilityAggressiveSolver : public VisibilitySolver<P, S>
    {
    public:
        VisibilityAggressiveSolver(VisibilityExactQuery_<P, S>* aSolver,
            S tolerance,
            bool detectApertureOnly);

        VisibilityResult resolve();
    private:
        void extractStabbingLines(PluckerPolyhedron<P>* myPolyhedron, PluckerPolytope<P>* aPolytope);

        bool mDetectApertureOnly;
        S mTolerance;
    };

    template<class P, class S>
    VisibilityAggressiveSolver<P, S>::VisibilityAggressiveSolver(VisibilityExactQuery_<P, S>* mQuery, S tolerance, bool detectApertureOnly)
        : VisibilitySolver<P, S>(mQuery),
        mTolerance(tolerance),
        mDetectApertureOnly(detectApertureOnly)
    {
    }


    template <class P, class S>
    VisibilityResult VisibilityAggressiveSolver<P, S>::resolve()
    {
        VisibilityResult myGlobalResult = HIDDEN;
        GeometryConvexPolygon* q0 = VisibilitySolver<P, S>::mQuery->getQueryPolygon(0);
        GeometryConvexPolygon* q1 = VisibilitySolver<P, S>::mQuery->getQueryPolygon(1);
  
        
     //   std::random_device rd;  // Will be used to obtain a seed for the random number engine
       // std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

        //std::uniform_real_distribution<> distribution(0.0, 1.0);

        std::vector<double> triangleFanArea0, triangleFanArea1;
        double polygonArea0 = MathGeometry::getTriangleFanAreas(*q0, triangleFanArea0);
        double polygonArea1 = MathGeometry::getTriangleFanAreas(*q1, triangleFanArea1);
        
        const size_t cumulativeProbabilitySize = 2048;
        std::vector<int> cumulativeProbability0, cumulativeProbability1;
        MathGeometry::computeCumulativeProbabilityLookupTable(cumulativeProbabilitySize,triangleFanArea0,cumulativeProbability0);
        MathGeometry::computeCumulativeProbabilityLookupTable(cumulativeProbabilitySize,triangleFanArea1,cumulativeProbability1);
        MathVector3d myBegin;
        MathVector3d myEnd;
        for (int i = 0; i < 2000; i++)
        {
            int triangle0  = 0;
            int triangle1  = 0;
        
            if (q0->getVertexCount() > 2)
            {
                int randomTriangleChoice0 = (int)(MathArithmetic<float>::getRandom() * (float)cumulativeProbabilitySize);
                triangle0 = cumulativeProbability0[randomTriangleChoice0];
            }
            if (q1->getVertexCount() > 2)
            {            
                int randomTriangleChoice1 = (int)(MathArithmetic<float>::getRandom() * (float)cumulativeProbabilitySize);
                 triangle1 = cumulativeProbability1[randomTriangleChoice1];
            }
         //   MathVector2d u1(distribution(gen), distribution(gen));
            {
                const std::vector<MathVector3d>& v0 = q0->getVertices();
                MathVector2d u0(MathArithmetic<double>::getRandom(),MathArithmetic<double>::getRandom());
                MathVector3d b0 = MathGeometry::sampleUniformTriangle(u0);                
                
                myBegin =      MathVector3d(b0.x * v0[0          ].x , b0.x * v0[0          ].y, b0.x * v0[0         ].z );
                if (v0.size() > 0)
                    myBegin += MathVector3d(b0.y * v0[triangle0+1].x , b0.y * v0[triangle0+1].y, b0.y * v0[triangle0+1].z );
                if (v0.size() > 1)
                    myBegin += MathVector3d(b0.z * v0[triangle0+2].x , b0.z * v0[triangle0+2].y, b0.z * v0[triangle0+2].z );
            }
            
            {
              //  MathVector2d u2(distribution(gen), distribution(gen));
                const std::vector<MathVector3d>& v1 = q1->getVertices();
  
                MathVector2d u1(MathArithmetic<double>::getRandom(),MathArithmetic<double>::getRandom());
                MathVector3d b1 = MathGeometry::sampleUniformTriangle(u1);
                myEnd = MathVector3d(b1.x * v1[0].x , b1.x * v1[0].y, b1.x * v1[0].z );
                if (v1.size() > 0)
                    myEnd += MathVector3d(b1.y * v1[triangle1+1].x , b1.y * v1[triangle1+1].y, b1.y * v1[triangle1+1].z );
                if (v1.size() > 1)
                    myEnd += MathVector3d(b1.z * v1[triangle1+2].x , b1.z * v1[triangle1+2].y, b1.z * v1[triangle1+2].z );
            }                  
               // bool hasIntersection = VisibilitySolver<P, S>::mQuery->hasSceneIntersection(myBegin,myEnd);
                bool hasIntersection = false;
                if (!hasIntersection)
                {                

                    if (VisibilitySolver<P, S>::mDebugger != nullptr)
                    {
                        VisibilitySolver<P, S>::mDebugger->addStabbingLine(convert<MathVector3f>(myBegin), convert<MathVector3f>(myEnd));
                    }

                    if (mDetectApertureOnly) 
                        return VISIBLE;
                    else
                        myGlobalResult = VISIBLE;   
                }
            
        }

        return myGlobalResult;

    }



}
