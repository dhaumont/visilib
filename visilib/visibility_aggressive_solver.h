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
#include "math_matrix_4.h"

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

    class GeometryConvexPolygonRandomSampler
    {
        public:
             GeometryConvexPolygonRandomSampler(const GeometryConvexPolygon& aPolygon,  size_t aCumulativeProbabilityTableSize = 2048)
             {
                initialize(aPolygon, aCumulativeProbabilityTableSize);
             }
             
             MathVector3d getSpatialSample()
             {
                int triangle  = 0;
                int vertexCount = mPolygon->getVertexCount();
                if (vertexCount > 2)
                {
                    int randomTriangleChoice = (int)(MathArithmetic<float>::getRandom() * (mCumulativeProbabilityTableSize-1));
                    triangle = mCumulativeProbabilityTable[randomTriangleChoice];
                }
                const std::vector<MathVector3d>& v = mPolygon->getVertices();
                MathVector2d u(getRandom(),getRandom());
                MathVector3d b = MathGeometry::uniformSampleTriangle(u);                
                
                MathVector3d myResult(b.x * v[0          ].x , b.x * v[0          ].y, b.x * v[0         ].z );
                if (vertexCount > 0)
                {   
                    myResult += MathVector3d(b.y * v[triangle+1].x , b.y * v[triangle+1].y, b.y * v[triangle+1].z );
                    if (vertexCount > 1)
                    {
                        myResult += MathVector3d(b.z * v[triangle+2].x , b.z * v[triangle+2].y, b.z * v[triangle+2].z );
                    }
                }
                return myResult;
             }


            MathVector3d getDirectionSample()
            {
                MathVector2d u(getRandom(),getRandom());
                MathVector3d d = MathGeometry::cosineSampleHemisphere(u);
                d = mLocalToWorldMatrix.multiply(d);        
                return d;
            }
            const double getPolygonArea() const
            {
                return mPolygonArea;
            }

            static size_t getSampleCount(double minimumHoleArea, double samplingArea, double confidenceValue, double errorMargin)
            {
                double p = 1 - minimumHoleArea / samplingArea;
                std::cout << "p:" << p << "; minimumHoleArea:" << minimumHoleArea << "; samplingArea:" <<  samplingArea << std::endl;
                //Normal distribution
                //double t = MathGeometry::getZValueForConfidenceLevel(confidenceValue, 0.0001);
                size_t count = 1;
                double cumulatedProbability = 1;
                while (cumulatedProbability >= (1.0 - confidenceValue))
                {
                    //std::cout << cumulatedProbability << " " << count << std::endl;
                    cumulatedProbability *= p;
                    count = count + 1;                    
                }
                std::cout << cumulatedProbability << " " << count << std::endl;
                //Cochran formulae
                //std::cout << "t: " << t << "; t * t * p * (1-p): " << t * t * p * (1-p) << "; errorMargin*errorMargin:" <<  errorMargin*errorMargin << std::endl;
                //size_t count = size_t(t * t * p * (1-p) / (errorMargin * errorMargin));
                
                std::cout << "COUNT: " << count << std::endl;
                
                return count;
            }
        private:

             double getRandom()
             {
                return MathArithmetic<double>::getRandom();
             }
            void initialize(const GeometryConvexPolygon& aPolygon, size_t aCumulativeProbabilityTableSize)
            {
                mCumulativeProbabilityTableSize = (float)aCumulativeProbabilityTableSize;
                mPolygon = &aPolygon;
                std::vector<double> myTriangleFanArea;
                mPolygonArea = MathGeometry::getTriangleFanAreas(aPolygon, myTriangleFanArea);
                MathGeometry::computeCumulativeProbabilityLookupTable(aCumulativeProbabilityTableSize, myTriangleFanArea, mCumulativeProbabilityTable);
                V_ASSERT(mCumulativeProbabilityTable.size() == aCumulativeProbabilityTableSize);

                MathPlane3d myPlane = aPolygon.getPlane();
                MathVector3d myNormal = myPlane.getNormal();
                myNormal.normalize();
                MathVector3d u,v;
                MathGeometry::getTangentBasis(myNormal, u, v);
                mLocalToWorldMatrix.setOrthogonalBasisInverseTranspose(u, v, myNormal);
                std::cout << mLocalToWorldMatrix;

                
                MathVector3d myZ(0,0,1);
                myZ = mLocalToWorldMatrix.multiply(myZ);                
            }
            std::vector<int> mCumulativeProbabilityTable;
            const GeometryConvexPolygon* mPolygon;
            float mCumulativeProbabilityTableSize;
            MathMatrixd mLocalToWorldMatrix;
            double mPolygonArea;
    };


    template <class P, class S>
    VisibilityResult VisibilityAggressiveSolver<P, S>::resolve()
    {
        VisibilityResult myGlobalResult = HIDDEN;
        const GeometryConvexPolygon& q0 = *VisibilitySolver<P, S>::mQuery->getQueryPolygon(0);
        const GeometryConvexPolygon& q1 = *VisibilitySolver<P, S>::mQuery->getQueryPolygon(1);
  
        GeometryConvexPolygonRandomSampler sampler0(q0, 2048);
        GeometryConvexPolygonRandomSampler sampler1(q1, 2048);
        
    
        double totalArea = sampler0.getPolygonArea() + sampler1.getPolygonArea();
        double polygonAreaRatio = sampler0.getPolygonArea() / (totalArea);
        size_t sampleCount = GeometryConvexPolygonRandomSampler::getSampleCount(0.0005, totalArea, 0.95, 0.05);
        
        for (size_t i = 0; i < sampleCount;)
        {
            bool isTriangle0WhichIsSampled = MathArithmetic<float>::getRandom() < polygonAreaRatio;
            
            GeometryConvexPolygonRandomSampler& sampler = isTriangle0WhichIsSampled ? sampler0 : sampler1;
            const GeometryConvexPolygon& myOppositePolyhon = isTriangle0WhichIsSampled ? q1 : q0;

            MathVector3d myBegin  = sampler.getSpatialSample();    
            MathVector3d myDirection = sampler.getDirectionSample();
            
            MathVector3d myEnd;  
            bool isSampleValid = false;
            
            if (MathGeometry::getPlaneIntersection(myOppositePolyhon.getPlane(), myBegin, myDirection, myEnd, MathArithmetic<double>::Tolerance()))
            {
                if (MathGeometry::isPointInsidePolygon(myOppositePolyhon,myEnd, MathArithmetic<double>::Tolerance()))
                {
                    isSampleValid = true;
                }
            }
                
            if (isSampleValid)
            {
                i+=1;
                bool hasIntersection = VisibilitySolver<P, S>::mQuery->hasSceneIntersection(myBegin,myEnd);
                
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
        }

        return myGlobalResult;

    }



}
