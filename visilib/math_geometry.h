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

#include "math_plucker_2.h"
#include "math_plucker_6.h"
#include "math_predicates.h"
#include "silhouette_mesh_face.h"
#include "geometry_position_type.h"
#include "geometry_convex_polygon.h"
#include "geometry_aabbox.h"
#include "geometry_ray.h"

namespace visilib
{
    class GeometryRay;
    class GeometryAABB;
    template<class S>
    class PluckerPolytope;
    template<class S>
    class PluckerPolyhedron;

    /** @brief Provides geometrical functions.*/

    class MathGeometry
    {
    public:

        /** @brief  Compute if a point is inside a convex polygon */
        static bool isPointInsidePolygon(const GeometryConvexPolygon& aPolygon, const MathVector3d& aPoint, double tolerance);

        static MathVector3d uniformSampleTriangle(const MathVector2d& u);
        static MathVector3d cosineSampleHemisphere(const MathVector2d& u);

        static void getTangentBasis(const MathVector3d& aUnitVector, MathVector3d& u, MathVector3d& v );
        /** @brief  Compute if a ray has an intersection with a triangle
        
        The ray triangle intersection is computed using code from "Watertight ray/triangle intersection" in Journal of Graphics Tools
        http://jcgt.org/published/0002/01/05/paper.pdf
        */
        template <class S> static bool hitsTriangle(const GeometryRay& ray, const MathVector3_<S>& v0, const MathVector3_<S>& v1, const MathVector3_<S>& v2);

        template <class S> static bool hitsCylinder(const GeometryRay& ray, const S& distance, const MathVector3_<S>& v0, const MathVector3_<S>& v1, const MathVector3_<S>& v2);

        /** @brief  Compute if a ray has an intersection with a bounding box 
        
        Plane/axis aligned box intersection code from David Eberly
        http://www.geometrictools.com/LibMathematics/Intersection/Intersection.html
        */
        static bool hitsBoundingBox(const GeometryRay& aRay, const GeometryAABB& boundingBox);

        /** @brief Compute the ray parameters for watertight ray-triangle intersection */
        static void computeRayParameters(const GeometryRay& aRay, double& Sx, double& Sy, double& Sz, int& kx, int& ky, int& kz);

        /** @brief  Test if an edge is potentially inside the shaft */
        static bool isEdgePotentiallyInsideShaft(const std::vector<MathPlane3d>& aShaftPlanes, const MathVector3f& a, const MathVector3f& b, bool exclude);

        /** @brief  Interpolate two vertices given two weights parameters 
        
        The method implement the method described in "Splitting aComplex of Convex Polytopes In Any Dimension" by Bajaj and Patucci
        http://www.cs.utexas.edu/~bajaj/cvc/papers/1996/conference/p88-bajaj.pdf
       
        @param anOffset1: the weight of the first vertex
        @param anOffset2: the weight of the second vertex
        @param aVertice1: the first vertex
        @param aVertice2: the second vertex
        @result aResult: the interpolated vertex computed by aResult = alpha * Vertice1 + beta * aVertice2
        */
        template<class P, class S> static P interpolate(S anOffset1, S anOffset2, const P& aVertice1, const P& aVertice2, S tolerance);

        /** @brief  Get the 3D anchor point of a Plucker point 
        
         p105 Jiry Bittner's PHD Thesis
         A (ax, ay, az) = Dl x Ll / ||Dl||^2 where Dl is the directional part and Ll the locational part
        */
        template<class P> static MathVector3d getAnchorPoint(const P& line);

        /** @brief Compute the intersections of an edge in Plucker space with the Plucker Quadric  
        @param v1: first vertex of the edge
        @param v2: first vertex of the edge
        @param p1: position of v1 with respect to the Plucker Quadric (ie the sign of v1.dot(v1))
        @param p2: position of v2 with respect to the Plucker Quadric (ie the sign of v2.dot(v2))
        @param result: the set of intersection points (0, 1 or 2) with the Plucker Quadric
        @param newtonRaphson: use an iterative Newton Raphson optimisation to increase the precision of the results
        */
        template<class P, class S> static bool findPluckerEdgeWithQuadricIntersection(const P& v1, const P& v2, GeometryPositionType p1, GeometryPositionType p2, std::vector<P>& result, bool newtonRaphson, S tolerance);

        /** @brief Compute the intersection of an edge in Plucker space with the Plucker Quadric  */
        template<class P, class S> static bool findPluckerEdgeWithQuadricIntersection(const P& v1, const P& v2, GeometryPositionType p1, GeometryPositionType p2, std::vector<P>& result, S tolerance);

        /** @brief Test if an edge has an intersection with the Plucker Quadric  */
        template<class P, class S> static bool hasPluckerEdgeWithQuadricIntersection(const P& v1, const P& v2, GeometryPositionType p1, GeometryPositionType p2, S tolerance);

        /** @brief Compute the supporting plane of a convex polygon */
        static MathPlane3d computePlane(const GeometryConvexPolygon& polygon);
     
        /** @brief Compute the supporting plane of a convex polygon using an approximate normal heuristic*/
        static MathPlane3d computePlaneFromApproximateNormal(const GeometryConvexPolygon& polygon, const MathVector3d& approximateNormal);

        /** @brief Compute the supporting plane of a triangle */
        template<class S>
        static MathPlane3_<S> computePlane(const MathVector3_<S>& v0, const MathVector3_<S>& v1, const MathVector3_<S>& v2);
    
        /** @brief  Clip a polygon following a plane equation, with an epsilon guard band */
        static bool clipWithGardBand(GeometryConvexPolygon& polygon, const MathPlane3d& aPlane, double anEpsilon);

        /** @brief Compute the 3D line of a Plucker point using two input 3D planes
        
        The method intersects the Plucker line with two planes, and return the two intersection points defining the 3D line.
        */
        template<class P> static std::pair<MathVector3d, MathVector3d> getBackTo3D(const P& aPoint, const MathPlane3d& aPlane1, const MathPlane3d& aPlane2);

        /** @brief Compute the 3D line of a Plucker point by computing intersection with a bounding box*/
        template<class P> static std::pair<MathVector3d, MathVector3d> getBackTo3D(const P& aPoint);

        /** @brief  Compute the intersection of a plane and a ray in 3D

         plane: ax+by+cz+d = 0
         line: P = aBegin + t * aDirection;
         intersection point: t = - (a * aBegin.x + b * aBegin.y + c * aBegin.z + d)/(a  * aDirection.x + b* aDirection.y + z * aDirection.z);
                               = - (myPlaneDirection.dotProduct(myDir) + w) / (myPlaneDirection.dotProduct(myDir));
        */
       template<class S>
        static bool getPlaneIntersection(const MathPlane3_<S>& plane, const MathVector3_<S>& aBegin, const MathVector3_<S>& aDirection, MathVector3_<S>& anIntersection, double epsilon);

        /** @brief Compute if a polygon is intersected by a plane */
        static bool intersect(const GeometryConvexPolygon& polygon, const MathPlane3d& aPlane);
       
        /** @brief  Compute the gravity center of a polygon */
        static MathVector3d getGravityCenter(const GeometryConvexPolygon& polygon);

        template<class S> static MathVector3_<S> getGravityCenter(const MathVector3_<S>& v0, const MathVector3_<S>& v1, const MathVector3_<S>& v2);

        static bool hasVertexOnFrontSide(const MathPlane3d& plane, SilhouetteMeshFace* face);

        /** @brief   */
        static bool hasVertexOutsidePlane(const MathPlane3d& plane, SilhouetteMeshFace* face);

        /** @brief  Compute a Plucker line that is representative of a polytope
        
        This procedure is not formally correct from a mathematical point of view: the gravity center in Pluker space do not correspond to a real line in 3D since the gravity center
         does not belongs to the Plucker quadric. However, it works well in practice to find an average polytope representative line
        */
        template<class P, class S> static P computeRepresentativeLine(PluckerPolytope<P>* polytope, PluckerPolyhedron<P>* polyhedron, S tolerance);

        /** @brief  Compute the intersection of an edge with a plane */
        static MathVector3_<double> getPlaneIntersectionWithEdge(const MathVector3_<double>& myV1, const MathVector3_<double>& myV2, const MathPlane3_<double>& aPlane);
   

        template<class P, class S> static P getProjectionOnQuadric(const P& line);
        template<class P, class S> static P getClosestQuadricPoint(const P& line);

        template<class P, class S>
        bool static isEdgeInsidePolytope(const MathVector3d& a, const MathVector3d& b, PluckerPolytope<P>* aPolytope, const MathVector3d& approximateNormal, PluckerPolyhedron<P>* polyhedron, S tolerance);
        
        template<class S>
        static bool isBoxInsideConvexHull(const MathVector3_<S>& AABBMin, const MathVector3_<S>& AABBMax, const std::vector<MathPlane3_<S> >& convexHullPlanes);

        static MathVector3d cartesianToSpherical(const MathVector3d& cartesian, bool computeLength = true);
        static MathVector3d sphericalToCartesian(const MathVector3d& spherical);

        static double getTriangleFanAreas(const GeometryConvexPolygon& aPolygon, std::vector<double>& triangleFanAreas);
        static void computeCumulativeProbabilityLookupTable(size_t lookupSize, const std::vector<double>& distribution, std::vector<int>& loopkupTable);
 };

    inline bool MathGeometry::isPointInsidePolygon(const GeometryConvexPolygon& aPolygon, const MathVector3d& aPoint, double tolerance)
    {
        if (aPolygon.getVertexCount() == 0)
            return false;

        if (aPolygon.getVertexCount() == 1)
        {
            MathVector3d r = aPoint - aPolygon.getVertex(0);
            return r.getSquaredNorm() < tolerance;
        }
        if (aPolygon.getVertexCount() == 2)
        {
            MathVector3d ab = aPolygon.getVertex(0) - aPolygon.getVertex(1);
            MathVector3d ac = aPolygon.getVertex(0) - aPoint;
            MathVector3d bc = aPolygon.getVertex(1) - aPoint;

            double test = MathArithmetic<double>::getAbs(ac.getNorm() + bc.getNorm() - ab.getNorm());
            //std::cout << "Test 2D: " << test << std::endl;
            return test < 1e-6; //TODO
        }
        if (MathArithmetic<double>::getAbs(aPolygon.getPlane().dot(aPoint)) > tolerance)
            return false;

        for (int i = 0; i < aPolygon.getVertexCount(); i++)
        {
            MathVector2i myEdge = aPolygon.getEdge(i);
            MathVector3d myNormal = MathVector3d::cross(aPolygon.getPlane().mNormal, aPolygon.getVertex(myEdge.x) - aPolygon.getVertex(myEdge.y));
            MathPlane3d myPlane(myNormal.x, myNormal.y, myNormal.z, -myNormal.dot(aPolygon.getVertex(i)));
            if (myPlane.dot(aPoint) > tolerance)
            {
                return false;
            }
        }

        return true;
    }

    inline bool MathGeometry::isEdgePotentiallyInsideShaft(const std::vector<MathPlane3d> & aShaftPlanes, const MathVector3f & a, const MathVector3f & b, bool exclude)
    {
        V_ASSERT(aShaftPlanes.size() > 0);

        //Approximate using a bounding sphere

        MathVector3f myCenter = a;
        myCenter += b;
        myCenter *= 0.5;

        MathVector3f radiusVector = b;
        radiusVector -= a;
        double radius = radiusVector.getNorm()*0.5;
        
        MathVector3d  myCenterD = convert<MathVector3d>(myCenter);
        for (size_t i = 0; i < aShaftPlanes.size(); i++)
        {
            double d = aShaftPlanes[i].dot(myCenterD);
    
            if (exclude)
            {
                if (d < -radius)
                {
                    return false;
                }
            }
            else
            {
                if (d <= -radius)
                {
                    return false;
                }
            }
        }
        return true;
    }

    inline bool MathGeometry::hasVertexOnFrontSide(const MathPlane3d & plane, SilhouetteMeshFace * face)
    {
        for (size_t i = 0; i < face->getVertexCount(); i++)
        {
            if (plane.dot(convert<MathVector3d>(face->getVertex(i))) > 0.0)//TODO ADD A TOLERANCE?
            {
                return true;
            }
        }

        return false;
    }

    inline MathVector3d MathGeometry::getGravityCenter(const GeometryConvexPolygon & polygon)
    {
        MathVector3d center;

        for (int i = 0; i < polygon.getVertexCount(); i++)
        {
            center += polygon.getVertex(i);
        }

        center *= (1. / (double)(polygon.getVertexCount()));

        return center;
    }


    template<class S>
    inline MathVector3_<S> MathGeometry::getGravityCenter(const MathVector3_<S>& v0, const MathVector3_<S>& v1, const MathVector3_<S> & v2)
    {
        MathVector3_<S> center(v0);        
        center+=v1;
        center+=v2;
        
        center *= (1. / 3.);

        return center;
    }

    template<class P>
    inline MathVector3d MathGeometry::getAnchorPoint(const P & line)
    {
        MathVector3d direction = convert<MathVector3d>(line.getDirection());
        MathVector3d location = convert<MathVector3d>(line.getLocation());
        MathVector3d result = MathVector3_<double>::cross(direction, location);
        double norm = direction.getSquaredNorm();

        if (MathArithmetic<double>::getAbs(norm - 1.0) >= MathArithmetic<double>::Tolerance())
        {
            if (norm <= 0)
            {
                std::cerr << "Can not normalize vector: " << std::endl;// << result.x << "," << result.y << "," << result.z << std::endl;
                result = MathVector3d::Zero();
                V_ASSERT(0);
            }
            else
            {
                result /= norm;
            }
        }

        return result;
    }

    template<class P, class S>
    inline P  MathGeometry::interpolate(S anOffset1, S anOffset2, const P & aVertice1, const P & aVertice2, S tolerance)
    {
        P  myVector1 = aVertice1;
        P  myVector2 = aVertice2;

        S lambdaV1 = MathArithmetic<S>::getAbs(anOffset1);
        S lambdaV2 = MathArithmetic<S>::getAbs(anOffset2);

        S eta = lambdaV2 - lambdaV1;
        S gamma;

        S sum;
        S sum2;

        if (eta < -tolerance)
        {
            eta = -eta;
            sum = lambdaV2;
            gamma = lambdaV2;
            sum2 = lambdaV2;
            sum += eta;
        }

        else
        {
            sum = lambdaV1;
            gamma = lambdaV1;
            sum2 = lambdaV1;
            myVector1 = aVertice2;
            sum += eta;
            myVector2 = aVertice1;
        }

        S alpha, beta;

        sum2 += sum;

        if (eta >= tolerance || gamma >= tolerance)
        {
            V_ASSERT(MathArithmetic<S>::getAbs(sum2) > tolerance);
            S div = 1.f / (sum2);
            alpha = gamma * div;
            beta = sum * div;
        }

        else
        {
            alpha = beta = 0.5f;
        }

        myVector1 *= alpha;
        myVector2 *= beta;
        myVector1 += myVector2;

        return myVector1;
    }

    template<>
    inline bool MathGeometry::findPluckerEdgeWithQuadricIntersection<MathPlucker6<float>, float>(const MathPlucker6<float> & v1, const MathPlucker6<float> & v2, GeometryPositionType p1, GeometryPositionType p2, std::vector<MathPlucker6<float> > & result, float tolerance)
    {
        return findPluckerEdgeWithQuadricIntersection<MathPlucker6<float>, float>(v1, v2, p1, p2, result, true, tolerance);
    }

    template<>
    inline bool MathGeometry::findPluckerEdgeWithQuadricIntersection<MathPlucker6<double>, double>(const MathPlucker6<double> & v1, const MathPlucker6<double> & v2, GeometryPositionType p1, GeometryPositionType p2, std::vector<MathPlucker6<double> > & result, double tolerance)
    {
        return findPluckerEdgeWithQuadricIntersection<MathPlucker6<double>, double>(v1, v2, p1, p2, result, true, tolerance);
    }

    template<>
    inline bool MathGeometry::findPluckerEdgeWithQuadricIntersection<MathPlucker2<double>, double>(const MathPlucker2<double> & v1, const MathPlucker2<double> & v2, GeometryPositionType p1, GeometryPositionType p2, std::vector<MathPlucker2<double> > & result, double tolerance)
    {
        return findPluckerEdgeWithQuadricIntersection<MathPlucker2<double>, double>(v1, v2, p1, p2, result, true, tolerance);
    }

#ifdef EXACT_ARITHMETIC
    template<>
    inline bool MathGeometry::findPluckerEdgeWithQuadricIntersection<MathPlucker6<exact>, exact>(const MathPlucker6<exact> & v1, const MathPlucker6<exact> & v2, GeometryPositionType p1, GeometryPositionType p2, std::vector<MathPlucker6<exact> > & result, exact tolerance)
    {
        return findPluckerEdgeWithQuadricIntersection<MathPlucker6<exact>, exact>(v1, v2, p1, p2, result, false, tolerance);
    }

#endif

    template<class P, class S>
    inline bool MathGeometry::findPluckerEdgeWithQuadricIntersection(const P & v1, const P & v2, GeometryPositionType p1, GeometryPositionType p2, std::vector<P> & result, bool newtonRaphson, S tolerance)
    {
        result.clear();
        
        if (p1 == ON_BOUNDARY)
        {
            result.push_back(v1);
        }
        if (p2 == ON_BOUNDARY)
        {
            result.push_back(v2);
        }
        if (result.size() > 0)
        {
            return true;
        }
        P  F = v2 - v1;

        S a = F.dot(F);
        S b = F.dot(v1);
        S c = v1.dot(v1);

        std::vector<S> ts;

        if (a <= tolerance && a >= -tolerance)
        {//a==0
            if (b <= tolerance && b >= -tolerance)
            {//b==0
                return false;
            }

            ts.push_back(-c / (2 * b));
        }
        else
        {
            S b_2 = b * b;
            if (a * a + b_2 + c * c <= tolerance)
            {
                return false;
            }
            S delta = b_2 - a * c;
            if (delta <= -tolerance)
            {
                return false; // delta < 0: no root
            }
            if (delta <= tolerance)
            {
                // delta is 0, single root
                ts.push_back(-b / a);
            }
            else
            {
                //general case: delta is positive, two roots
                S sqrt_delta = MathArithmetic<S>::getSqrt(delta);
                V_ASSERT(MathArithmetic<S>::getAbs(a) >= tolerance);

                ts.push_back((-b + sqrt_delta) / a);
                ts.push_back((-b - sqrt_delta) / a);
            }
        }
        for (size_t i = 0; i < ts.size(); i++)
        {
            S t = ts[i];

            if (t >= -tolerance && t <= 1 + tolerance)
            {
                if (newtonRaphson)
                {
                    S at, b2;
                    b2 = 2 * b;
                    S fxn = (a * t * t + b2 * t + c);

                    int iter = 0;

                    // Newton - Raphson iteration scheme : xn + 1 = xn - f(xn) / f'(xn)
                    //f (t)	 = a t*t + 2 b t + c
                    //f'(t)  = 2 at + 2b
                    const double toleranceLocal = 1e-18;
                    while (MathArithmetic<S>::getAbs(fxn) >= toleranceLocal && iter < 20)
                    {
                        at = a * t;
                        fxn = (at * t + b2 * t + c);
                        S fdxn = (2 * at + b2);
                        t = t - fxn / fdxn;
                        iter++;
                    }
                }
                if (t >= -tolerance && t <= 1 + tolerance)
                {
                    result.push_back((F * t + v1));
                }
            }
        }

#if _DEBUG
        for (size_t i = 0; i < result.size(); i++)
        {
            V_ASSERT(result[i].dot(result[i]) < tolerance);
        }
#endif
        return result.size() > 0;
    }

    template<class P, class S>
    inline bool MathGeometry::hasPluckerEdgeWithQuadricIntersection(const P & v1, const P & v2, GeometryPositionType p1, GeometryPositionType p2, S tolerance)
    {
        if (p1 == ON_BOUNDARY || p2 == ON_BOUNDARY)
        {
            return true;
        }

        P  F = v2 - v1;

        S a = F.dot(F);
        S b = F.dot(v1);
        S c = v1.dot(v1);

        std::vector<S> ts;

        if (a <= tolerance && a >= -tolerance)
        {
            if (b <= tolerance && b >= -tolerance)
            {
                return false;
            }

            S t = (-c / (2 * b));
            return (t >= -tolerance && t <= 1 + tolerance);
        }
        else
        {
            S b_2 = b * b;
            if (a * a + b_2 + c * c <= tolerance)
            {
                return false;
            }
            S delta = b_2 - a * c;
            if (delta <= -tolerance)
            {
                return false; // delta < 0: no root
            }
            if (delta <= tolerance)
            {
                // delta is 0, single root
                S t = (-b / a);
                return (t >= -tolerance && t <= 1 + tolerance);
            }
            else
            {
                //general case: delta is positive, two roots
                S sqrt_delta = MathArithmetic<S>::getSqrt(delta);
                V_ASSERT(MathArithmetic<S>::getAbs(a) >= tolerance);

                S t = (-b + sqrt_delta) / a;
                if (t >= -tolerance && t <= 1 + tolerance)
                    return true;
                t = (-b - sqrt_delta) / a;
                if (t >= -tolerance && t <= 1 + tolerance)
                    return true;

                return false;
            }
        }

        return false;
    }

    /** @brief
    //Code from http://www.terathon.com/code/clipping.html
    // A little bit modified: we translate the plane of a distance anEpsilon
    // (to remove a supplementary strip band)
    */

    inline bool MathGeometry::clipWithGardBand(GeometryConvexPolygon & aPolygon, const MathPlane3d & aPlane, double anEpsilon)
    {
        MathPlane3d myPlane(aPlane);
        myPlane.d -= anEpsilon;

        enum
        {
            polygonInterior = 1,
            polygonBoundary = 0,
            polygonExterior = -1
        };

        std::vector<signed char> location(aPolygon.getVertexCount());

        int positive = 0;
        int negative = 0;

        for (size_t a = 0; a < aPolygon.getVertexCount(); a++)
        {
            double d = myPlane.dot(aPolygon.getVertex(a));
            if (d > MathArithmetic<double>::Tolerance())
            {
                location[a] = polygonInterior;
                positive++;
            }
            else if (d < -MathArithmetic<double>::Tolerance())
            {
                location[a] = polygonExterior;
                negative++;
            }
            else
            {
                location[a] = polygonBoundary;
            }
        }

        if (negative == 0)
        {
            return true;
        }
        else if (positive == 0)
        {
            aPolygon.clear();
            return false;
        }
        //	std::cout<<"CLIP -the query polygon"<<std::endl;
        size_t previous = aPolygon.getVertexCount() - 1;

        std::vector<MathVector3d> result;

        for (int index = 0; index < aPolygon.getVertexCount(); index++)
        {
            long loc = location[index];
            if (loc == polygonExterior)
            {
                if (location[previous] == polygonInterior)
                {
                    const MathVector3d& v1 = aPolygon.getVertex(previous);
                    const MathVector3d& v2 = aPolygon.getVertex(index);
                    MathVector3d myRes = getPlaneIntersectionWithEdge(v1, v2, myPlane);
                    result.push_back(myRes);
                }
            }
            else
            {
                const MathVector3d& v1 = aPolygon.getVertex(index);
                if ((loc == polygonInterior) && (location[previous] == polygonExterior))
                {
                    const MathVector3d& v2 = aPolygon.getVertex(previous);
                    MathVector3d myRes = getPlaneIntersectionWithEdge(v1, v2, myPlane);
                    result.push_back(myRes);
                }

                result.push_back(v1);
            }

            previous = index;
        }
        //TODO TEST! (not tested)
        aPolygon = GeometryConvexPolygon(result);

        V_ASSERT(aPolygon.getVertexCount() > 1);
        return (aPolygon.getVertexCount() > 1);
    }

    inline MathPlane3d MathGeometry::computePlane(const GeometryConvexPolygon & polygon)
    {
        bool edgeNotFound = true;
        size_t i = 1;
        MathVector3d myNormal;

        while (edgeNotFound && i < polygon.getVertexCount() - 1)
        {
            myNormal = MathVector3d::cross(polygon.getVertex(i) - polygon.getVertex(i - 1), polygon.getVertex(i + 1) - polygon.getVertex(i));

            edgeNotFound = (myNormal.getSquaredNorm() < MathArithmetic<double>::Tolerance());
            i++;
        }

        MathVector3d axes[3] = { MathVector3d::X(), MathVector3d::Y(), MathVector3d::Z() };

        if (edgeNotFound && polygon.getVertexCount() > 1)
        {
            MathVector3d myEdge1 = polygon.getVertex(0) - polygon.getVertex(1);

            while (edgeNotFound && i < 3)
            {
                myNormal = MathVector3d::cross(myEdge1, axes[i]);
                edgeNotFound = (myNormal.getSquaredNorm() < MathArithmetic<double>::Tolerance());
                i++;
            }
        }
        if (edgeNotFound)
        {
            myNormal = axes[0];
            edgeNotFound = false;
        }
        V_ASSERT(!edgeNotFound);
        myNormal.normalize();

        return MathPlane3d(myNormal.x, myNormal.y, myNormal.z, -myNormal.dot(polygon.getVertex(0)));
    }

    template<class S>
    inline MathPlane3_<S> MathGeometry::computePlane(const MathVector3_<S>& v0, const MathVector3_<S>& v1, const MathVector3_<S>& v2)
    {
        MathVector3_<S> myNormal = MathVector3_<S>::cross(v1 - v0, v2 - v0);
        myNormal.normalize();

        return MathPlane3_<S>(myNormal.x, myNormal.y, myNormal.z, -myNormal.dot(v0));
    }

    inline MathPlane3d MathGeometry::computePlaneFromApproximateNormal(const GeometryConvexPolygon & aPolygon, const MathVector3d & approximateNormal)
    {
        V_ASSERT(aPolygon.getVertexCount() >= 2);

        MathVector3d edge = aPolygon.getVertex(1) - aPolygon.getVertex(0);

        MathVector3d e1 = MathVector3d::cross(approximateNormal, edge);

        MathVector3d myNormal = MathVector3d::cross(edge, e1);
        myNormal.normalize();

        return MathPlane3d(myNormal.x, myNormal.y, myNormal.z, -myNormal.dot(aPolygon.getVertex(0)));
    }

    template<class P>
    inline std::pair<MathVector3d, MathVector3d> MathGeometry::getBackTo3D(const P & aPoint, const MathPlane3d & aPlane1, const MathPlane3d & aPlane2)
    {
        MathVector3d myBegin = MathGeometry::getAnchorPoint(aPoint);
        MathVector3d myDir = convert<MathVector3d>(aPoint.getDirection());
        V_ASSERT(!myDir.isZero(MathArithmetic<double>::Tolerance()));
        double myMax = myDir.normalize();
        V_ASSERT(myMax > MathArithmetic<double>::Tolerance());
        MathVector3d myResult1;
        MathVector3d myResult2;
        bool result1 = getPlaneIntersection(aPlane1, myBegin, myDir, myResult1, MathArithmetic<double>::Tolerance());
        bool result2 = getPlaneIntersection(aPlane2, myBegin, myDir, myResult2, MathArithmetic<double>::Tolerance());
        V_ASSERT(result1 && result2);

        return std::pair<MathVector3d, MathVector3d>(myResult1, myResult2);
    }

    template<class P>
    inline std::pair<MathVector3d, MathVector3d> MathGeometry::getBackTo3D(const P & aPoint)
    {
        MathVector3d axis[3] = { MathVector3d::X(), MathVector3d::Y(), MathVector3d::Z() };
        int	condition = 0;
        MathVector3d myDir = convert<MathVector3d>(aPoint.getDirection());

        double max = MathArithmetic<double>::getAbs(myDir.dot(axis[0]));

        for (int i = 1; i < 3; i++)
        {
            double dx = MathArithmetic<double>::getAbs(myDir.dot(axis[i]));
            if (dx > max)
            {
                max = dx;
                condition = i;
            }
        }

        MathPlane3d plane1(axis[condition].x, axis[condition].y, axis[condition].z, 1);
        MathPlane3d plane2(axis[condition].x, axis[condition].y, axis[condition].z, -1);

        return getBackTo3D(aPoint, plane1, plane2);
    }

    template<class S>
    inline bool MathGeometry::getPlaneIntersection(const MathPlane3_<S> & plane, const MathVector3_<S> & aBegin, const MathVector3_<S> & aDirection, MathVector3_<S> & anIntersection, double tolerance)
    {
        anIntersection = aBegin;
        MathVector3_<S> myDir(aDirection);

        S lambda = -(plane.dot(aBegin));
        S div = plane.mNormal.dot(myDir);

        if (MathArithmetic<S>::getAbs(div) > tolerance)
        {
            lambda /= div;
            myDir *= lambda;
            anIntersection += myDir;
            V_ASSERT(MathArithmetic<S>::getAbs(plane.dot(anIntersection)) <= tolerance);

            return true;
        }

        return false;
    }

    template<class S>
    inline bool MathGeometry::hitsTriangle(const GeometryRay & ray, const MathVector3_<S> & v0, const MathVector3_<S> & v1, const MathVector3_<S> & v2)
    {
        MathVector3_<S> A = v0 - convert<MathVector3_<S> >(ray.getStart());
        MathVector3_<S> B = v1 - convert<MathVector3_<S> >(ray.getStart());
        MathVector3_<S> C = v2 - convert<MathVector3_<S> >(ray.getStart());

        S Sx;	S Sy;	S Sz;
        int kx;	int ky;	int kz;
        ray.get(Sx, Sy, Sz, kx, ky, kz);

        const S Ax = A[kx] - Sx * A[kz];
        const S Ay = A[ky] - Sy * A[kz];
        const S Bx = B[kx] - Sx * B[kz];
        const S By = B[ky] - Sy * B[kz];
        const S Cx = C[kx] - Sx * C[kz];
        const S Cy = C[ky] - Sy * C[kz];

        S U = Cx * By - Cy * Bx;
        S V = Ax * Cy - Ay * Cx;
        S W = Bx * Ay - By * Ax;

        // Fall back to test against edges using S precision.
        if (U == 0.0 || V == 0.0 || W == 0.0) {
            S CxBy = static_cast<S>(Cx) * static_cast<S>(By);
            S CyBx = static_cast<S>(Cy) * static_cast<S>(Bx);
            U = static_cast<S>(CxBy - CyBx);

            S AxCy = static_cast<S>(Ax) * static_cast<S>(Cy);
            S AyCx = static_cast<S>(Ay) * static_cast<S>(Cx);
            V = static_cast<S>(AxCy - AyCx);

            S BxAy = static_cast<S>(Bx) * static_cast<S>(Ay);
            S ByAx = static_cast<S>(By) * static_cast<S>(Ax);
            W = static_cast<S>(BxAy - ByAx);
        }

        /*if (trace_options_.cull_back_face) {
        if (U < 0.0 || V < 0.0 || W < 0.0) return false;
        }
        */
        if ((U < 0.0 || V < 0.0 || W < 0.0) && (U > 0.0 || V > 0.0 || W > 0.0)) {
            return false;
        }

        S det = U + V + W;
        if (det == 0.0) return false;

        /*
        const S Az = Sz * A[kz];
        const S Bz = Sz * B[kz];
        const S Cz = Sz * C[kz];
        const S D = U * Az + V * Bz + W * Cz;

        const S rcpDet = S(1.0) / det;

        S tt = D * rcpDet;
        if (tt > (*t_inout))
        {
        return false;
        }
        */
        return true;
    }


    template<class S>
    inline bool MathGeometry::hitsCylinder(const GeometryRay & ray, const S& distance, 
                                           const MathVector3_<S> & v0,
                                           const MathVector3_<S> & v1, 
                                           const MathVector3_<S> & v2)
    {
        const S epsilon = MathArithmetic<S>::Tolerance();
        MathPlane3_<S> myPlane;
        myPlane = MathPlane3_<S>::computeTrianglePlane(v0,v1,v2, epsilon);

        MathVector3_<S> myRayPiercingPointInPlane;
        bool hit = getPlaneIntersection(myPlane, ray.getStart(), ray.getDirection(), myRayPiercingPointInPlane, epsilon);
        if (!hit)
            return false;
        
        MathVector3_<S> gravityCenter = getGravityCenter(v0,v1,v2);
        
        S enclosingCircleRadius = (gravityCenter-v0).getSquaredNorm();
        enclosingCircleRadius = std::max(enclosingCircleRadius, (gravityCenter -v1).getSquaredNorm());
        enclosingCircleRadius = std::max(enclosingCircleRadius, (gravityCenter -v2).getSquaredNorm());
        enclosingCircleRadius = MathArithmetic<S>::getSqrt(enclosingCircleRadius);

        S d = distance + enclosingCircleRadius;
        return (myRayPiercingPointInPlane - gravityCenter).getSquaredNorm() < d * d;
    }

    template<class P, class S>
    bool  MathGeometry::isEdgeInsidePolytope(const MathVector3d& a, const MathVector3d& b, PluckerPolytope<P>* aPolytope, const MathVector3d& approximateNormal, PluckerPolyhedron<P>* polyhedron, S tolerance)
    {
        MathVector3d e = b;
        e -= a;
        MathVector3d e1 = MathVector3d::cross(approximateNormal, e);
        e1.normalize();
        double d = e.getNorm();
        e1 *= d;
        // Try to split the polytope by the hyperplane of the candidate edge
        bool hasPointInside1 = false;
        bool hasPointInside2 = false;

        {
            P  myHyperplane(a, a + e1);

            for (auto v : aPolytope->getVertices())
            {
                GeometryPositionType position = MathPredicates::getVertexPlaneRelativePosition(myHyperplane,polyhedron->get(v), tolerance);
                if (position == ON_NEGATIVE_SIDE || position == ON_BOUNDARY)
                {
                    hasPointInside1 = true;
                    break;
                }

            }
        }

        {
            P  myHyperplane(b, b + e1);

            for (auto v : aPolytope->getVertices())
            {
                GeometryPositionType position = MathPredicates::getVertexPlaneRelativePosition(myHyperplane, polyhedron->get(v), tolerance);
                if (position == ON_POSITIVE_SIDE || position == ON_BOUNDARY)
                {
                    hasPointInside2 = true;
                    break;
                }
            }
        }

        return hasPointInside1 && hasPointInside2;
    }

    template<class P, class S>
    inline P MathGeometry::computeRepresentativeLine(PluckerPolytope<P> * polytope, PluckerPolyhedron<P> * polyhedron, S tolerance)
    {   
        P  myGravityCenterImaginary = P::Zero();
        auto myVertices = polytope->getVertices();
        
        for (auto iter = myVertices.begin(); iter != myVertices.end(); iter++)
        {
            size_t v = *iter;
            myGravityCenterImaginary += polyhedron->get(v);
        }
        
        P myGravityCenterReal = getProjectionOnQuadric<P,S>(myGravityCenterImaginary);
//        P myGravityCenterReal = getClosestQuadricPoint<P,S>(myGravityCenterImaginary);
        V_ASSERT(MathPredicates::getQuadricRelativePosition(myGravityCenterReal, tolerance) == ON_BOUNDARY);
        return myGravityCenterReal;
    }

    inline void MathGeometry::computeRayParameters(const GeometryRay & aRay, double& Sx, double& Sy, double& Sz, int& kx, int& ky, int& kz)
    {
        const MathVector3f& myDirection = aRay.getDirection();

        kz = 0;
        double absDir = fabs(myDirection[0]);
        if (absDir < fabs(myDirection[1]))
        {
            kz = 1;
            absDir = fabs(myDirection[1]);
        }
        if (absDir < fabs(myDirection[2]))
        {
            kz = 2;
            absDir = fabs(myDirection[2]);
        }

        kx = kz + 1;
        if (kx == 3) kx = 0;
        ky = kx + 1;
        if (ky == 3) ky = 0;

        // Swap kx and ky dimention to preserve widing direction of triangles.
        if (myDirection[kz] < 0.0f) std::swap(kx, ky);

        // Claculate shear constants.
        Sx = myDirection[kx] / myDirection[kz];
        Sy = myDirection[ky] / myDirection[kz];
        Sz = 1.0f / myDirection[kz];
    }

    inline bool MathGeometry::hitsBoundingBox(const GeometryRay & aRay, const GeometryAABB & boundingBox)
    {
        MathVector3f myBoxHalfSize = (boundingBox.getMax() - boundingBox.getMin()) * 0.5f;
        MathVector3f myBoxCenter = boundingBox.getMin() + myBoxHalfSize;

        float fWdU[3], fAWdU[3], fDdU[3], fADdU[3], fAWxDdU[3], fRhs;

        MathVector3f kDiff = aRay.getStart() - myBoxCenter;

        fWdU[0] = aRay.getDirection().x;
        fAWdU[0] = std::fabs(fWdU[0]);
        fDdU[0] = kDiff.x;
        fADdU[0] = std::fabs(fDdU[0]);
        if (fADdU[0] > myBoxHalfSize.x && fDdU[0] * fWdU[0] >= 0.0f)
            return false;

        fWdU[1] = aRay.getDirection().y;
        fAWdU[1] = std::fabs(fWdU[1]);
        fDdU[1] = kDiff.y;
        fADdU[1] = std::fabs(fDdU[1]);
        if (fADdU[1] > myBoxHalfSize.y && fDdU[1] * fWdU[1] >= 0.0f)
            return false;

        fWdU[2] = aRay.getDirection().z;
        fAWdU[2] = std::fabs(fWdU[2]);
        fDdU[2] = kDiff.z;
        fADdU[2] = std::fabs(fDdU[2]);
        if (fADdU[2] > myBoxHalfSize.z && fDdU[2] * fWdU[2] >= 0.0f)
            return false;

        MathVector3f kWxD = MathVector3f::cross(aRay.getDirection(), kDiff);

        fAWxDdU[0] = std::fabs(kWxD.x);
        fRhs = myBoxHalfSize.y * fAWdU[2] + myBoxHalfSize.z * fAWdU[1];
        if (fAWxDdU[0] > fRhs)
            return false;

        fAWxDdU[1] = std::fabs(kWxD.y);
        fRhs = myBoxHalfSize.x * fAWdU[2] + myBoxHalfSize.z * fAWdU[0];
        if (fAWxDdU[1] > fRhs)
            return false;

        fAWxDdU[2] = std::fabs(kWxD.z);
        fRhs = myBoxHalfSize.x * fAWdU[1] + myBoxHalfSize.y * fAWdU[0];
        if (fAWxDdU[2] > fRhs)
            return false;

        return true;
    }

    inline MathVector3_<double> MathGeometry::getPlaneIntersectionWithEdge(const MathVector3_<double> & myV1, const MathVector3_<double> & myV2, const MathPlane3_<double> & aPlane)
    {
        double offset1 = aPlane.dot(myV1);
        double offset2 = aPlane.dot(myV2);

        return interpolate<MathVector3d, double>(offset1, offset2, myV1, myV2, MathArithmetic<double>::Tolerance());
    }

    inline bool MathGeometry::hasVertexOutsidePlane(const MathPlane3d & plane, SilhouetteMeshFace * face)
    {
        for (int i = 0; i < face->getVertexCount(); i++)
        {
            if (std::fabs(plane.dot(convert<MathVector3d>(face->getVertex(i)))) >= MathArithmetic<double>::Tolerance())
            {
                return true;
            }
        }

        return false;
    }

    template<class P, class S> 
    inline P MathGeometry::getProjectionOnQuadric(const P& line)
    {
        const MathVector3_<S>& a = line.getDirection();
        const MathVector3_<S>& b = line.getLocation();
        S p = a.dot(b);
        S q = a.getSquaredNorm() + b.getSquaredNorm();

        S mu = 2 * p / ( q+ MathArithmetic<S>::getSqrt(q*q-4*p*p));

        return P((a.x - mu * b.x), (a.y - mu * b.y), (a.z - mu * b.z), (b.x - mu * a.x), (b.y - mu * a.y), (b.z - mu * a.z));
    }

    /** @brief Return the closest Plucker point on the Plucker quadric 
        Solve the Plucker correction problem as described in the paper
        "Plücker Correction Problem: Analysis and Improvements in Efficiency"*/        
    template<class P, class S> 
    inline P MathGeometry::getClosestQuadricPoint(const P& line)
    {            
        const MathVector3_<S>& a = line.getDirection();
        const MathVector3_<S>& b = line.getLocation();
        const S a1 = a.x, a2 = a.y, a3 = a.z;
        const S b1 = b.x,  b2 = b.y, b3 = b.z;

        const S p = a1*b1 + a2*b2 + a3*b3;
        const S q = a1*a1 + a2*a2 + a3*a3 + b1*b1 + b2*b2 + b3*b3;
        const S mu = 2*p/(q+MathArithmetic<S>::getSqrt(q*q-4*p*p));
        const S u_ = 1/(1-mu*mu);
        
        const S x1 = (a1-mu*b1)*u_, x2 = (a2-mu*b2)*u_, x3 = (a3-mu*b3)*u_;
        const S y1 = (b1-mu*a1)*u_, y2 = (b2-mu*a2)*u_, y3 = (b3-mu*a3)*u_;
        
        return P(x1, x2, x3, y1, y2, y3);
    }

    template<class S>
    inline bool isBoxInsideConvexHull(const MathVector3_<S>& AABBMin, const MathVector3_<S>& AABBMax, const std::vector<MathPlane3_<S> >& convexHullPlanes)
    {
            for (const auto& plane: convexHullPlanes)
            {
                double d = std::max(AABBMin.x * plane.getNormal().x, AABBMax.x * plane.getNormal().x) 
                         + std::max(AABBMin.y * plane.getNormal().y, AABBMax.y * plane.getNormal().y) 
                         + std::max(AABBMin.z * plane.getNormal().z, AABBMax.z * plane.getNormal().z) 
                         + plane.d;

                if (d>0)  return false;
            }
            return true;
    }


    inline void GeometryConvexPolygon::computePlane()
    {
        mPlane = MathGeometry::computePlane(*this);
    }

    inline void GeometryRay::initialize(const MathVector3f & aStart, const MathVector3f & aDirection)
    {
        mStart = aStart;
        mDirection = aDirection;

        MathGeometry::computeRayParameters(*this, Sx, Sy, Sz, kx, ky, kz);
    }

    template<>
    inline bool MathGeometry::getPlaneIntersection<float>(const MathPlane3f & plane, const MathVector3f & aBegin, const MathVector3f & aDirection, MathVector3f & anIntersection, double tolerance)
    {
        MathVector3d myIntersection;
        bool result = MathGeometry::getPlaneIntersection<double>(convert<MathPlane3d>(plane),
                                                                 convert<MathVector3d>(aBegin),
                                                                 convert<MathVector3d>(aDirection),
                                                                 myIntersection,
                                                                tolerance);
        anIntersection = convert<MathVector3f>(myIntersection);
        return result;
    }

    inline MathVector3d MathGeometry::uniformSampleTriangle(const MathVector2d& u) 
    {
        double b0, b1;
        if (u.x < u.y) 
        {
            b0 = u.x / 2;
            b1 = u.y - b0;
        } 
        else 
        {
            b1 = u.y / 2;
            b0 = u.x - b1;
        }
        return MathVector3d(b0, b1, 1 - b0 - b1);
    }

    inline MathVector3d MathGeometry::cosineSampleHemisphere(const MathVector2d& u)
    {
        const double r = sqrt(u.x);
        const double theta = 2 * M_PI * u.y;
 
        const double x = r * cos(theta);
        const double y = r * sin(theta);
 
        return MathVector3d(x, y, sqrt(std::max(0.0, 1. - u.x)));
    }
    inline MathVector3d MathGeometry::cartesianToSpherical(const MathVector3d& cartesian, bool computeLength)
    {
        double theta = atan2(cartesian.y, cartesian. x);
        double phi = acos(cartesian.z);
        double r = 1;
        if (computeLength)
        {
            r = cartesian.getNorm();
        }

        return MathVector3d(r,theta,phi);
    }
    
    inline MathVector3d MathGeometry::sphericalToCartesian(const MathVector3d& spherical)
    {
        double sinPhi = sin(spherical.z);
        double cosPhi = cos(spherical.z);
        double sinTheta = sin(spherical.y);
        double cosTheta = cos(spherical.y);
        double r = spherical.z;

        
        return MathVector3d(r * cosTheta * sinPhi, r * sinTheta * sinPhi, r * cosPhi);
    }

    inline double MathGeometry::getTriangleFanAreas(const GeometryConvexPolygon& aPolygon, std::vector<double>& triangleFanAreas)
    {
        double totalArea = 0.0;    
        const MathVector3d& v0 = aPolygon.getVertex(0);
        
        MathVector3d myNormal;
        for (size_t i = 1; i < aPolygon.getVertexCount()-1;i++)
        {
            MathVector3d myEdge0(aPolygon.getVertex(i),v0);
            MathVector3d myEdge1(aPolygon.getVertex(i+1),v0);
            myNormal = MathVector3d::cross(myEdge0,myEdge1);
            double triangleArea = myNormal.getNorm();
            totalArea += triangleArea;
            triangleFanAreas.push_back(triangleArea);
        }

        return totalArea;
    }

    inline void MathGeometry::computeCumulativeProbabilityLookupTable(size_t lookupSize, const std::vector<double>& distribution, std::vector<int>& loopkupTable)
    {
        if (distribution.empty())
        {
            return;
        }
        double distributionSum = 0;
        
        for (auto p : distribution)
        {
            distributionSum += p;
        }
        const double increment = 1.0 / (double)lookupSize;
        double currentIndex = 0;
        double cumulatedProbability = 0;
        double cumulatedProbabilityEndOfInterval = distribution[0];
        loopkupTable.resize(lookupSize);
        
        for (size_t i = 0; i < lookupSize; i++)
        {
            if (cumulatedProbability > cumulatedProbabilityEndOfInterval)
            {
                currentIndex++;
                cumulatedProbabilityEndOfInterval += distribution[currentIndex];
            }
           loopkupTable[i] = currentIndex;
           cumulatedProbability += (increment * distributionSum);         
        }
    }

    inline void MathGeometry::getTangentBasis(const MathVector3d& aUnitVector, MathVector3d& u, MathVector3d& v )
    {
    //    V_ASSERT(MathArithmetic<double>::getAbs(aUnitVector.getSquaredNorm() - 1.0) < MathArithmetic<double>::Tolerance );
        MathVector3d av(MathArithmetic<double>::getAbs(aUnitVector.x),MathArithmetic<double>::getAbs(aUnitVector.y),MathArithmetic<double>::getAbs(aUnitVector.z) );
        
        MathVector3d mainAxis;
        if (av.x >= av.z && av.x >= av.y)
        {
            mainAxis.x = 0; mainAxis.y = 1;mainAxis.z = 0;
        }
        else if (av.y >= av.x && av.y >= av.z)
        {
            mainAxis.x = 0; mainAxis.y = 0;mainAxis.z = 1;
        }
        else
        {   
            V_ASSERT(av.z >= av.x && av.z >= av.x);
            
            mainAxis.x = 1; mainAxis.y = 0;mainAxis.z = 0;
        }
        
        u = MathVector3d::cross(aUnitVector, mainAxis);
        u.normalize();

        v = MathVector3d::cross(aUnitVector,u);
        v.normalize();        

        std::cout << "n:" << aUnitVector <<  "main:" << mainAxis  << "; u: " << u << "; v: " << v << std::endl; 
        

    }
}
