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
#include "visibility_aperture_finder.h"

#include "geometry_convex_polygon.h"
#include "silhouette_mesh_face.h"
#include "math_plucker_2.h"
#include "math_predicates.h"
#include "plucker_polytope.h"
#include "plucker_polytope_splitter.h"
#include "visibility_exact_query.h"
#include "silhouette.h"
#include <fstream>
#include <sstream>
#include <map>
namespace visilib
{
    template<class P>
    class  PluckerPolytope;

    /** @brief Exact visibility determination algorithm based on an aperture finder heuristic. The algorithm stops as soon as an aperture was detected.

    This class manage the CSG operations sheduling in Plucker space.
    It implements an algorithm that searches if at least on aperture exist between the two source polygons, such that at least one line that is intersecting
    the two source polygons is not blocked by any occluder. */

    template<class P, class S>
    class VisibilityApertureFinder : public VisibilitySolver<P, S>
    {
    public:
        VisibilityApertureFinder(VisibilityExactQuery_<P, S>* aSolver,
            bool normalization,
            S tolerance,
            bool detectApertureOnly);

        VisibilityResult resolve();
    private:
        VisibilityResult resolveInternal(VisibilityResult& aGlobalResult, PluckerPolytope<P,S>* aPolytope, const std::string& occlusionTreeNodeSymbol, const std::vector<Silhouette<P>*>& anOccluders, const std::vector<P>& aPolytopeLines, int depth);
        void resize(size_t myInitiaLineCount, PluckerPolytope<P,S>* aPolytope);
        void extractStabbingLines(PluckerPolyhedron<P>* myPolyhedron, PluckerPolytope<P,S>* aPolytope);
            

        bool mNormalization;
        bool mDetectApertureOnly;
        S mTolerance;
    };

    template<class P, class S>
    VisibilityApertureFinder<P, S>::VisibilityApertureFinder(VisibilityExactQuery_<P, S>* mQuery, bool normalization, S tolerance, bool detectApertureOnly)
        : VisibilitySolver<P, S>(mQuery),
        mNormalization(normalization),
        mTolerance(tolerance),
        mDetectApertureOnly(detectApertureOnly)
    {
    }

    template <class P, class S>
    VisibilityResult VisibilityApertureFinder<P, S>::resolve()
    {
        VisibilityResult myGlobalResult = HIDDEN;

        resolveInternal(myGlobalResult, reinterpret_cast<PluckerPolytope<P,S>*>(VisibilitySolver<P, S>::mQuery->getComplex()->getRoot()), "*", std::vector<Silhouette<P>*>(), std::vector<P>(), 0);

        return myGlobalResult;
    }

    template<class P, class S>
    VisibilityResult VisibilityApertureFinder<P, S>::resolveInternal(VisibilityResult& aGlobalResult, PluckerPolytope<P,S>* aPolytope, const std::string& occlusionTreeNodeSymbol, const std::vector<Silhouette<P>*>& anOccluders, const std::vector<P>& aPolytopeLines, int aDepth)
    {
        size_t myInitiaLineCount = myPolyhedron->getLinesCount();

#ifdef OUTPUT_DEBUG_FILE
        std::ofstream& debugOutput = VisibilitySolver<P, S>::mDebugger->getDebugOutput();    
        V_LOG(debugOutput, "Resolve internal ", occlusionTreeNodeSymbol);
        aPolytope->outputProperties(debugOutput, myPolyhedron);
#endif
        if (aDepth > 2000)
        {
            std::cerr << "Overflow detected: more than 200 recursive calls...." << std::endl;
#ifdef OUTPUT_DEBUG_FILE
            V_LOG(debugOutput, "RESULT FAILURE: Overflow detected : more than 2000 recursive calls....", occlusionTreeNodeSymbol);
#endif
            aGlobalResult = FAILURE;
            return FAILURE;
        }
        {
            HelperScopedTimer timer(VisibilitySolver<P, S>::mQuery->getStatistic(), STABBING_LINE_EXTRACTION);
            aPolytope->computeEdgesIntersectingQuadric(mTolerance);
        }
        if (!aPolytope->containsRealLines())
        {
#ifdef OUTPUT_DEBUG_FILE
            V_LOG(debugOutput, "RESULT HIDDEN: polytope does not contains real line", occlusionTreeNodeSymbol);
#endif
            return HIDDEN;
        }

        V_ASSERT(aPolytope->isValid(myPolyhedron, mNormalization, mTolerance));

        bool hasRay = false;
        std::vector<Silhouette<P>*> myOccluders = anOccluders;
        std::vector<P> polytopeLines = aPolytopeLines;
        if (myOccluders.empty())
        {
            if (!VisibilitySolver<P, S>::mQuery->collectAllOccluders(aPolytope, myPolyhedron, myOccluders, polytopeLines))
            {
                // Early stop - an aperture has been found: the source polygons are mutually visible
#ifdef OUTPUT_DEBUG_FILE
                V_LOG(debugOutput, "RAY: IS VISIBLE", occlusionTreeNodeSymbol);
#endif
                hasRay = true;
                aGlobalResult = VISIBLE;
                if (mDetectApertureOnly)
                {
                    return VISIBLE;
                }
            }
            else
            {
#ifdef OUTPUT_DEBUG_FILE
                V_LOG(debugOutput, "RAY: IS OCCLUDED", occlusionTreeNodeSymbol);
#endif
            }
        }

{
        HelperScopedTimer timer(VisibilitySolver<P, S>::mQuery->getStatistic(), OCCLUDER_TREATMENT);
        if (VisibilitySolver<P, S>::mQuery->isOccluded(aPolytope, myOccluders, polytopeLines))
        {
#ifdef OUTPUT_DEBUG_FILE
            V_LOG(debugOutput, "COUNT: IS OCCLUDED", occlusionTreeNodeSymbol);
#endif
            return HIDDEN;
        }
        else
        {
#ifdef OUTPUT_DEBUG_FILE
            V_LOG(debugOutput, "COUNT: maybe VISIBLE", occlusionTreeNodeSymbol);
#endif

        }
}
        Silhouette<P>* mySilhouette<P> = nullptr;
        size_t mySilhouette<P>EdgeIndex = 0;

        bool hasEdge = false;
        bool isVisible = false;
        hasEdge = VisibilitySolver<P, S>::mQuery->findNextEdge(mySilhouette<P>EdgeIndex, mySilhouette<P>, aPolytope, occlusionTreeNodeSymbol);

        if (hasEdge) // a candidate edge has been found. we will split the polytope with the hyperplane of this edge.
        {
            Silhouette<P>Edge& myVisibilitySilhouette<P>Edge = mySilhouette<P>->getEdge(mySilhouette<P>EdgeIndex);
            SilhouetteMeshFace* face = myVisibilitySilhouette<P>Edge.mFace;

            V_ASSERT(myVisibilitySilhouette<P>Edge.mIsActive);
            // deactivate the candidate edge for further recursion
            mySilhouette<P>->setEdgeActive(mySilhouette<P>EdgeIndex, false);

            MathVector2i edge = face->getEdge(myVisibilitySilhouette<P>Edge.mEdgeIndex);

            MathVector3d a = convert<MathVector3d>(face->getVertex(edge.x));
            MathVector3d b = convert<MathVector3d>(face->getVertex(edge.y));

            bool intersect = false;

            {
                HelperScopedTimer timer(VisibilitySolver<P, S>::mQuery->getStatistic(), OCCLUDER_TREATMENT);
                intersect = MathGeometry::isEdgeInsidePolytope(a, b, aPolytope, VisibilitySolver<P, S>::mQuery->getApproximateNormal(),  mTolerance);
            }
            if (intersect)
            {
                if (VisibilitySolver<P, S>::mDebugger != nullptr)
                {
                    VisibilitySolver<P, S>::mDebugger->addRemovedEdge(face->getVertex(edge.x), face->getVertex(edge.y));
                }

                size_t myPolyhedronFace = myVisibilitySilhouette<P>Edge.mHyperPlaneIndex;

                if (myPolyhedronFace == 0)
                {
                    P myHyperplane(a, b);
                    if (mNormalization)
                    {
                        myHyperplane = myHyperplane.getNormalized();
                    }
                    myPolyhedronFace = myPolyhedron->add(myHyperplane, ON_BOUNDARY, mNormalization, mTolerance);
                    myVisibilitySilhouette<P>Edge.mHyperPlaneIndex = myPolyhedronFace;
                }

                P myHyperplane = myPolyhedron->get(myPolyhedronFace);
                
                GeometryPositionType myResult = ON_NEGATIVE_SIDE;

                PluckerPolytope<P,S>* myPolytopeLeft = new PluckerPolytope<P,S>();
                PluckerPolytope<P,S>* myPolytopeRight = new PluckerPolytope<P,S>();

                {
                    HelperScopedTimer timer(VisibilitySolver<P, S>::mQuery->getStatistic(), POLYTOPE_SPLIT);
                    VisibilitySolver<P, S>::mQuery->getStatistic()->inc(POLYTOPE_SPLIT_COUNT);

#ifdef OUTPUT_DEBUG_FILE
                    V_LOG(debugOutput, "PERFORM THE SPLIT", occlusionTreeNodeSymbol);
#endif
                    myResult = PluckerPolytopeSplitter<P, S>::split(myPolyhedron, myHyperplane, aPolytope, myPolytopeLeft, myPolytopeRight, myPolyhedronFace, mNormalization, mTolerance);

                    if (VisibilitySolver<P, S>::mQuery->getStatistic()->get(POLYTOPE_SPLIT_COUNT) % 10000 == 0)
                    {
                        VisibilitySolver<P, S>::mQuery->getStatistic()->displayCounts();
                        std::cout << std::endl;
                    }
                }

                std::vector <PluckerPolytope<P,S>*> myPolytopes;
                std::vector <std::string> postFix;
                std::vector<bool> reuseOccluders;

                if (myResult == ON_BOUNDARY)
                {
                    V_ASSERT(myPolytopeLeft && myPolytopeRight);
#ifdef OUTPUT_DEBUG_FILE
                    V_LOG(debugOutput, "SPLIT SUCCESS ->recurse on Left and right childs", occlusionTreeNodeSymbol);
#endif
                    auto position = MathPredicates::getVertexPlaneRelativePosition(myHyperplane, polytopeLines[0], mTolerance);

                    // left split
                    reuseOccluders.push_back(position != ON_POSITIVE_SIDE);
                    myPolytopes.push_back(myPolytopeLeft);
                    postFix.push_back("L");

                    // right split
                    reuseOccluders.push_back(position != ON_NEGATIVE_SIDE);
                    myPolytopes.push_back(myPolytopeRight);
                    postFix.push_back("R");
                }
                else
                {
#ifdef OUTPUT_DEBUG_FILE
                    V_LOG(debugOutput, "NO SPLIT OCCURS -> recurse again on the same polytope", occlusionTreeNodeSymbol);
#endif
                    reuseOccluders.push_back(true);
                    myPolytopes.push_back(aPolytope);
                    postFix.push_back("*");
                    delete myPolytopeLeft;  myPolytopeLeft = nullptr;
                    delete myPolytopeRight; myPolytopeRight = nullptr;
                }

                for (size_t i = 0; i < myPolytopes.size(); i++)
                {
                    bool hasBeenAdded = false;

                    if (myResult == ON_BOUNDARY && i == 0 || myResult == ON_NEGATIVE_SIDE)
                    {
                        mySilhouette<P>->pushEdgeProcessed(mySilhouette<P>EdgeIndex);
                        hasBeenAdded = true;
                    }


                    std::stringstream ss;
#ifdef OUTPUT_DEBUG_FILE
                    ss << occlusionTreeNodeSymbol << "|" << face->getFaceIndex() << "-" << myVisibilitySilhouette<P>Edge.mEdgeIndex << postFix[i];
#endif
                    VisibilityResult result;

                    if (reuseOccluders[i])
                    {
                        result = resolveInternal(aGlobalResult, myPolytopes[i], ss.str(), myOccluders, polytopeLines, aDepth + 1);
                    }
                    else
                    {
                        result = resolveInternal(aGlobalResult, myPolytopes[i], ss.str(), std::vector<Silhouette<P>*>(), std::vector<P>(), aDepth + 1);
                    }
                    if (result == FAILURE)
                    {
                        aGlobalResult = FAILURE;
                        return FAILURE;
                    }
                    if (result == VISIBLE)
                    {
                        aGlobalResult = VISIBLE;
                        if (mDetectApertureOnly) // Early stop - an aperture has been found
                        {
#ifdef OUTPUT_DEBUG_FILE
                            V_LOG(debugOutput, "EARLY STOP - aperture found");
#endif
                            delete myPolytopeLeft;  myPolytopeLeft = nullptr;
                            delete myPolytopeRight; myPolytopeRight = nullptr;

                            return VISIBLE;
                        }
                    }
                    if (hasBeenAdded)
                        mySilhouette<P>->popEdgeProcessed(mySilhouette<P>EdgeIndex);

                }
                delete myPolytopeLeft;  myPolytopeLeft = nullptr;
                delete myPolytopeRight; myPolytopeRight = nullptr;

            }
            else
            {

                std::stringstream ss;
#ifdef OUTPUT_DEBUG_FILE              
                ss << occlusionTreeNodeSymbol << "|" << face->getFaceIndex() << "-" << myVisibilitySilhouette<P>Edge.mEdgeIndex << "*";
#endif
                VisibilityResult result = resolveInternal(aGlobalResult, aPolytope, ss.str(), anOccluders, aPolytopeLines, aDepth + 1);
                if (result == FAILURE)
                    return result;

                if (result == VISIBLE)
                {
                    aGlobalResult = VISIBLE;
                    if (mDetectApertureOnly) // Early stop - an aperture has been found
                    {
                        return VISIBLE;
                    }
                }
            }
            // reactivate the  edge

            mySilhouette<P>->setEdgeActive(mySilhouette<P>EdgeIndex, true);
        }

        // No valid candidate edge has been found in the occluder set: all the set of lines represented by the polytope are blocked by at least one occluder: the recursion can stops.
        // When it is the case for all the sub-polytopes created during the recursive traversal, the two polygons are mutually hidden.

        if (!mDetectApertureOnly)
        { 
            if (!hasEdge)
            {
                //We found a final polytope visible
#ifdef OUTPUT_DEBUG_FILE
                V_LOG(debugOutput, "visible polytope found: exporting stabbing lines", occlusionTreeNodeSymbol);
#endif
                extractStabbingLines(myPolyhedron, aPolytope);
                aGlobalResult = VISIBLE;
                if (mDetectApertureOnly) // Early stop - an aperture has been found
                {
                    return VISIBLE;
                }
            }
        }

        if (myPolyhedron->getLinesCount() != myInitiaLineCount)
        {
       //     resize(myInitiaLineCount, myPolyhedron, aPolytope);
        }
        return UNKNOWN;

    }

    template<class P, class S>
    void VisibilityApertureFinder<P, S>::resize(size_t anInitiaLineCount, PluckerPolytope<P,S>* aPolytope)
    {
        HelperScopedTimer timer(VisibilitySolver<P, S>::mQuery->getStatistic(), POLYTOPE_SPLIT);

        myPolyhedron->resize(anInitiaLineCount);
        for (auto v : aPolytope->getVertices())
        {
            std::vector<size_t>& facets = myPolyhedron->getFacetsDescription(v);
            for (size_t i = 0; i < facets.size(); i++)
            {
                if (facets[i] >= anInitiaLineCount)
                {
                    facets.resize(i);
                }
            }
        }
    }

    template<class P, class S>
    void VisibilityApertureFinder<P, S>::extractStabbingLines(PluckerPolytope<P,S>* aPolytope)
    {
        HelperScopedTimer timer(VisibilitySolver<P, S>::mQuery->getStatistic(), STABBING_LINE_EXTRACTION);

        MathPlane3d aPlane0 = VisibilitySolver<P, S>::mQuery->getQueryPolygon(0)->getPlane();
        MathPlane3d aPlane1 = VisibilitySolver<P, S>::mQuery->getQueryPolygon(1)->getPlane();

        if (aPolytope->getExtremalStabbingLinesCount() == 0)
        {
            aPolytope->computeExtremalStabbingLines(myPolyhedron, mTolerance);
        }

        if (VisibilitySolver<P, S>::mDebugger != nullptr)
        {  
            std::vector<std::pair<MathVector3d, MathVector3d>> lines;
            aPolytope->getExtremalStabbingLinesBackTo3D(lines, aPlane0, aPlane1);
            for (auto line : lines)
            {
                VisibilitySolver<P, S>::mDebugger->addExtremalStabbingLine(convert<MathVector3f>(line.first), convert<MathVector3f>(line.second));
            }
        }

    }


}
