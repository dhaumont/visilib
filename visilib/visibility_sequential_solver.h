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

#include "geometry_convex_polygon.h"
#include "silhouette_mesh_face.h"
#include "math_plucker_2.h"
#include "math_predicates.h"
#include "plucker_polyhedron.h"
#include "plucker_polytope.h"
#include "plucker_polytope_complex.h"
#include "plucker_polytope_splitter.h"
#include "visibility_exact_query.h"
#include "silhouette.h"
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
namespace visilib
{
    template <class P>
    class PluckerPolytope;

    /** @brief  */

    template <class P, class S>
    class VisibilitySequentialSolver : public VisibilitySolver<P, S>
    {
    public:
        VisibilitySequentialSolver(VisibilityExactQuery_<P, S> *aSolver,
                                   bool normalization,
                                   S tolerance,
                                   bool detectApertureOnly);

        VisibilityResult resolve();

    private:
        PolytopePositionType classifyRelativeToSilhouette(PluckerPolytope<P> *myPolytopeToClassify,
                                          PluckerPolyhedron<P> *myPolyhedron,
                                          const std::vector<SilhouetteEdge> &edges
                                          );

        void splitWithAllSilhouetteEdges(std::vector<int> &workingArray,
                                         std::unordered_set<PluckerPolytope<P> *> &toDelete,
                                         const std::vector<SilhouetteEdge> &edges);

        void createMissingHyperplanes(PluckerPolyhedron<P> *aPolyhedron,
                                      std::vector<SilhouetteEdge> &edges);

        void classifyAllPolytopes(
            const std::vector<SilhouetteEdge> &edges,
            PluckerPolytopeComplex<P> *complex,
            PluckerPolyhedron<P> *myPolyhedron,
            std::vector<int> &intersectedPolytopes);

        void treatIntersectedPolytopes(
            const std::vector<SilhouetteEdge> &edges,
            PluckerPolytopeComplex<P> *complex,
            PluckerPolyhedron<P> *myPolyhedron,
            const std::vector<int> &intersectedPolytopes);

        bool mNormalization;
        bool mDetectApertureOnly;
        S mTolerance;
    };

    template <class P, class S>
    VisibilitySequentialSolver<P, S>::VisibilitySequentialSolver(VisibilityExactQuery_<P, S> *mQuery, bool normalization, S tolerance, bool detectApertureOnly)
        : VisibilitySolver<P, S>(mQuery),
          mNormalization(normalization),
          mTolerance(tolerance),
          mDetectApertureOnly(detectApertureOnly)
    {
    }

    template <class P, class S>
    void VisibilitySequentialSolver<P, S>::createMissingHyperplanes(PluckerPolyhedron<P> *aPolyhedron,
                                                                    std::vector<SilhouetteEdge> &edges)
    {
        for (SilhouetteEdge &myVisibilitySilhouetteEdge : edges)
        {
            size_t myPolyhedronFace = myVisibilitySilhouetteEdge.mHyperPlaneIndex;

            if (myPolyhedronFace == 0)
            {
                SilhouetteMeshFace *face = myVisibilitySilhouetteEdge.mFace;

                MathVector2i edge = face->getEdge(myVisibilitySilhouetteEdge.mEdgeIndex);

                MathVector3d a = convert<MathVector3d>(face->getVertex(edge.x));
                MathVector3d b = convert<MathVector3d>(face->getVertex(edge.y));

                P myHyperplane(a, b);
                if (mNormalization)
                {
                    myHyperplane = myHyperplane.getNormalized();
                }
                myPolyhedronFace = aPolyhedron->add(myHyperplane, ON_BOUNDARY, mNormalization, mTolerance);
                myVisibilitySilhouetteEdge.mHyperPlaneIndex = myPolyhedronFace;
            }
        }
    }

    template <class P, class S>
    PolytopePositionType VisibilitySequentialSolver<P, S>::classifyRelativeToSilhouette(PluckerPolytope<P> *aPolytope,
                                                                        PluckerPolyhedron<P> *aPolyhedron,
                                                                        const std::vector<SilhouetteEdge> &edges)
    {    
        bool intersect = false;
        bool positive = false;
        bool negative = false;
        bool boundary = false;
        
        if (!aPolytope->containsRealLines())
        {
            V_TRACE(" polytope  " << aPolytope << " is blocked ");
            V_TRACE("does not containsreallines: " << aPolytope);
            return BLOCKED;
        }
        if (edges.size() == 0)
        {
             V_TRACE(" polytope  " << aPolytope << " is blocked because no silhouette edge ");
            return BLOCKED;
        }
        V_TRACE("prepare to classify " << aPolytope << " against " << edges.size() <<" edges");
        for (auto myVisibilitySilhouetteEdge : edges)
        {            
            size_t myPolyhedronFace = myVisibilitySilhouetteEdge.mHyperPlaneIndex;
            V_ASSERT(myPolyhedronFace != 0);
            P myHyperplane = aPolyhedron->get(myPolyhedronFace);
            V_TRACE("Hyperplane: " << myVisibilitySilhouetteEdge.mHyperPlaneIndex << ", " << myHyperplane);
            auto myResult = MathPredicates::getRelativePosition(aPolytope, aPolyhedron, myHyperplane, mTolerance);
            V_TRACE(" classification of " << aPolytope << " against facet " << myPolyhedronFace << ":" <<myResult);
            switch (myResult)
            {
               case ON_BOUNDARY:       boundary = true;   break;
               case ON_BOTH_SIDES:     intersect = true;  break;
               case ON_NEGATIVE_SIDE:  negative = true;   break;
               case ON_POSITIVE_SIDE:  positive = true;   break;          
            }
            
        }

        if (!intersect)
        {
            if (positive && !negative)
            {
                return BLOCKED;
            }
            else if (negative && !positive)
            {
                return OUTSIDE;
            }
         }
        return INTERSECTS;       
    }

    template <class P, class S>
    void VisibilitySequentialSolver<P, S>::splitWithAllSilhouetteEdges(std::vector<int> &workingArray,
                                                                       std::unordered_set<PluckerPolytope<P> *> &toDelete,
                                                                       const std::vector<SilhouetteEdge> &edges)
    {
        PluckerPolytopeComplex<P> *complex = VisibilitySolver<P, S>::mQuery->getComplex();
        PluckerPolyhedron<P> *myPolyhedron = reinterpret_cast<PluckerPolyhedron<P> *>(VisibilitySolver<P, S>::mQuery->getComplex()->getPolyhedron());
         V_TRACE("Preparing to split polytopes with " << edges.size() << " silhouette edges");
         
        for (const SilhouetteEdge &myVisibilitySilhouetteEdge : edges)
        {
            SilhouetteMeshFace *face = myVisibilitySilhouetteEdge.mFace;

            V_ASSERT(myVisibilitySilhouetteEdge.mIsActive);

            MathVector2i edge = face->getEdge(myVisibilitySilhouetteEdge.mEdgeIndex);

            // Create the Plucker representation of the edge and add hyperplane of the edge and add it to the polyhedron

            size_t myPolyhedronFace = myVisibilitySilhouetteEdge.mHyperPlaneIndex;
            V_ASSERT(myPolyhedronFace != 0);

            if (VisibilitySolver<P, S>::mDebugger != nullptr)
            {
                VisibilitySolver<P, S>::mDebugger->addRemovedEdge(face->getVertex(edge.x), face->getVertex(edge.y));
            }

            P myHyperplane = myPolyhedron->get(myPolyhedronFace);
           V_TRACE("Preparing to split " << workingArray.size() << " polytopes");
                  
             
            GeometryPositionType myResult = ON_UNKNOWN_POSITION;
            for (int i = workingArray.size() - 1; i >= 0; i--)
            {
                PluckerPolytope<P> *myPolytopeToSplit = complex->getPolytope(workingArray[i]);
                V_TRACE("Split  " << myPolytopeToSplit);
         
                PluckerPolytope<P> *myPolytopeLeft = new PluckerPolytope<P>();
                PluckerPolytope<P> *myPolytopeRight = new PluckerPolytope<P>();
                V_TRACE("Creating myPolytopeLeft: " << myPolytopeLeft);
                V_TRACE("Creating myPolytopeRight: " << myPolytopeRight);
                {
                    HelperScopedTimer timer(VisibilitySolver<P, S>::mQuery->getStatistic(), POLYTOPE_SPLIT);
                    VisibilitySolver<P, S>::mQuery->getStatistic()->inc(POLYTOPE_SPLIT_COUNT);

                    myResult = PluckerPolytopeSplitter<P, S>::split(myPolyhedron, myHyperplane, myPolytopeToSplit, myPolytopeLeft, myPolytopeRight, myPolyhedronFace, mNormalization, mTolerance);

                    if (VisibilitySolver<P, S>::mQuery->getStatistic()->get(POLYTOPE_SPLIT_COUNT) % 10000 == 0)
                    {
                        VisibilitySolver<P, S>::mQuery->getStatistic()->displayCounts();
                        V_TRACE("");
                    }
                }

                  V_TRACE("result:  " << myResult);
         
                if (myResult == ON_BOTH_SIDES)
                {
                    V_ASSERT(MathPredicates::getRelativePosition(myPolytopeLeft, myPolyhedron, myHyperplane, mTolerance) == ON_NEGATIVE_SIDE);
                    V_ASSERT(MathPredicates::getRelativePosition(myPolytopeRight, myPolyhedron, myHyperplane, mTolerance) == ON_POSITIVE_SIDE);
                    size_t size_before = complex->getPolytopeCount();

                    myPolytopeLeft->computeEdgesIntersectingQuadric(myPolyhedron, mTolerance);
                    myPolytopeRight->computeEdgesIntersectingQuadric(myPolyhedron, mTolerance);
                 
                    if (myPolytopeLeft->containsRealLines())
                    {
                        complex->appendPolytope(myPolytopeLeft);
                    }
                    else
                    {
                        delete myPolytopeLeft;
                        myPolytopeLeft = nullptr;
                    }

                    if (myPolytopeRight->containsRealLines())
                    {
                        complex->appendPolytope(myPolytopeRight);
                    }
                    else
                    {
                        delete myPolytopeRight;
                        myPolytopeRight = nullptr;
                    }
                    toDelete.insert(myPolytopeToSplit);                    
                    V_TRACE("toDelete.insert:" << myPolytopeToSplit << " at index " << workingArray[i]);
                    
                    if (myPolytopeLeft && myPolytopeRight)
                    {
                        workingArray[i] = size_before;
                        workingArray.push_back(size_before + 1);
                        V_ASSERT(complex->getPolytope(workingArray[i]) == myPolytopeLeft);
                        V_ASSERT(complex->getPolytope(workingArray[workingArray.size() - 1]) == myPolytopeRight);
                    }
                    else if (myPolytopeLeft)
                    {
                        workingArray[i] = size_before;
                        V_ASSERT(complex->getPolytope(workingArray[i]) == myPolytopeLeft);
                      
                    }
                    else if (myPolytopeRight)
                    {
                       workingArray[i] = size_before;   
                        V_ASSERT(complex->getPolytope(workingArray[i]) == myPolytopeRight);
                    }
                    else
                    {
                        V_ASSERT(0);
                    }
                }
                else
                {
                    V_ASSERT(MathPredicates::getRelativePosition(myPolytopeToSplit, myPolyhedron, myHyperplane, mTolerance) == myResult);
                    delete myPolytopeLeft;
                    myPolytopeLeft = nullptr;
                    delete myPolytopeRight;
                    myPolytopeRight = nullptr;
                    V_TRACE("Delete left and right:" << myPolytopeLeft << " " << myPolytopeRight);
                  
                }
            }
        }
    }

    template <class P, class S>
    void VisibilitySequentialSolver<P, S>::classifyAllPolytopes(
        const std::vector<SilhouetteEdge> &edges,
        PluckerPolytopeComplex<P> *complex,
        PluckerPolyhedron<P> *myPolyhedron,
        std::vector<int> &intersectedPolytopes)
    {

          V_TRACE("Preparing to classify " << complex->getPolytopeCount() << " polytopes");
         
        for (int polytopeIndex = 0; polytopeIndex < complex->getPolytopeCount(); polytopeIndex++)
        {
            PluckerPolytope<P> *myPolytopeToClassify = complex->getPolytope(polytopeIndex);
            V_TRACE("classify " << myPolytopeToClassify << " at index " << polytopeIndex);
 
            PolytopePositionType position = classifyRelativeToSilhouette(myPolytopeToClassify, myPolyhedron, edges);

            if (position == INTERSECTS)
            {
                V_TRACE("adding to intersectedPolytopes "  << myPolytopeToClassify << " at " <<  polytopeIndex);
                intersectedPolytopes.push_back(polytopeIndex);
            }
            else if (position == BLOCKED)
            {
               V_TRACE("classifyRelativeToSilhouette is blocked: delete "  << myPolytopeToClassify << " at " <<  polytopeIndex);
               

                delete myPolytopeToClassify;
                complex->setPolytope(polytopeIndex, NULL);
            }
            else 
            {
                V_ASSERT(position == OUTSIDE);
                 V_TRACE("classifyRelativeToSilhouette is outside - skip: "  << myPolytopeToClassify << " at " <<  polytopeIndex);
             
            }
        }
    }

    template <class P, class S>
    void VisibilitySequentialSolver<P, S>::treatIntersectedPolytopes(
        const std::vector<SilhouetteEdge> &edges,
        PluckerPolytopeComplex<P> *complex,
        PluckerPolyhedron<P> *myPolyhedron,
        const std::vector<int> &intersectedPolytopes)
    {
        if (intersectedPolytopes.empty())
        {
               V_TRACE("treatIntersectedPolytopes : no intersectedpolytopes - done ");
         
        
            return;
        }
        V_TRACE("treatIntersectedPolytopes: prepare to treat " <<  intersectedPolytopes.size() << " intersected polytopes");
     
        for (int intersectedIndex = intersectedPolytopes.size() - 1; intersectedIndex >= 0; intersectedIndex--)
        {
            
            size_t originalPolytopeCount = complex->getPolytopeCount();
             V_TRACE("originalPolytopeCount stored as  " <<  originalPolytopeCount);
     
            std::unordered_set<PluckerPolytope<P> *> toDelete;
            std::vector<int> workingArray;
            int originalPolytopeIndex = intersectedPolytopes[intersectedIndex];

            workingArray.push_back(originalPolytopeIndex);


            V_TRACE("before splitWithAllSilhouetteEdgesm working array is:  " <<  workingArray.size());

            splitWithAllSilhouetteEdges(workingArray, toDelete, edges);

            bool atLeastOneNewPolytopeIsBlocked = false;
     
            V_TRACE("after splitWithAllSilhouetteEdges, working array is: " <<  workingArray.size());

            for (int i = 0; i < workingArray.size(); i++)
            {
                int myPolytopeToBlockIndex = workingArray[i];
                PluckerPolytope<P> *myPolytopeToBlock = complex->getPolytope(myPolytopeToBlockIndex);

                V_TRACE("computeEdgesIntersectingQuadric for polytope " <<  myPolytopeToBlock << " at index " << myPolytopeToBlockIndex);

                myPolytopeToBlock->computeEdgesIntersectingQuadric(myPolyhedron, mTolerance);
               
                 V_TRACE("classifyRelativeToSilhouette for polytope " <<  myPolytopeToBlock << " at index " << myPolytopeToBlockIndex);

                PolytopePositionType position = classifyRelativeToSilhouette(myPolytopeToBlock, myPolyhedron, edges);

                //V_ASSERT(position != INTERSECTS);
                if (position == BLOCKED)
                {

                 V_TRACE(myPolytopeToBlock <<  " isBlocked");

                    if (toDelete.find(myPolytopeToBlock) != toDelete.end())
                    {
                        V_TRACE("Double delete:" << myPolytopeToBlock << " " << myPolytopeToBlockIndex);
                    }
                     V_TRACE("toDelete.insert:" << myPolytopeToBlock << " at index " << workingArray[i]);
                  
                    toDelete.insert(myPolytopeToBlock);
                    atLeastOneNewPolytopeIsBlocked = true;
                    complex->setPolytope(myPolytopeToBlockIndex, NULL);
                    workingArray[i] = -1;
                }
            }
            
            if (atLeastOneNewPolytopeIsBlocked)
            {
                 V_TRACE("atLeastOneNewPolytopeIsBlocked");

               
                 V_TRACE("set polytope at index " << originalPolytopeIndex<< " is NULL ");

                 V_TRACE("Start iterating over polytopescount to set the deleted one to NULL ");

                for (int i=0; i < complex->getPolytopeCount(); i++)
                {
                    PluckerPolytope<P>* myPolytope = complex->getPolytope(i,false);
                    if (toDelete.find(myPolytope) != toDelete.end())
                    {
                        complex->setPolytope(i,NULL);
                    }
                }
               
                for (auto polytope: toDelete)
                {
                    V_TRACE("delete polytopes from delete: " << polytope);

                    delete polytope;
                }
            }
            else
            {
                V_TRACE("!atLeastOneNewPolytopeIsBlocked: revert change ");

                V_TRACE(" complex->getPolytopeCount() : " << complex->getPolytopeCount()  << " and  originalPolytopeCount:" << originalPolytopeCount);

                while (complex->getPolytopeCount() != originalPolytopeCount)
                {
                    V_TRACE("removelast");

                    complex->removeLast();
                     V_TRACE(" complex->getPolytopeCount() : " << complex->getPolytopeCount()  << " and  originalPolytopeCount:" << originalPolytopeCount);

         
                }
            }
        }
    }

    template <class P, class S>
    VisibilityResult VisibilitySequentialSolver<P, S>::resolve()
    {
        PluckerPolytopeComplex<P> *complex = VisibilitySolver<P, S>::mQuery->getComplex();
        PluckerPolyhedron<P> *myPolyhedron = reinterpret_cast<PluckerPolyhedron<P> *>(VisibilitySolver<P, S>::mQuery->getComplex()->getPolyhedron());

        const std::unordered_set<Silhouette *> &mySilhouettes = VisibilitySolver<P, S>::mQuery->getSilhouettes();

        for (auto iter = mySilhouettes.begin(); iter != mySilhouettes.end(); iter++)
        {
            Silhouette *mySilhouette = (*iter);
            auto &edges = mySilhouette->getEdges();
            createMissingHyperplanes(myPolyhedron, edges);

            std::vector<int> intersectedPolytopes;
            classifyAllPolytopes(edges, complex, myPolyhedron, intersectedPolytopes);

            treatIntersectedPolytopes(edges, complex, myPolyhedron, intersectedPolytopes);

            complex->repack();
        }
        V_TRACE("End of procedure:" << complex->getPolytopeCount() << " left");
        if (complex->getPolytopeCount() == 0)
        {
             V_TRACE("HIDDEN: done");
       
            return HIDDEN;
        }
         V_TRACE("visible: prepare to extract stabbinglines");
       
        for (int i = 0; i < complex->getPolytopeCount(); i++)
        {
            PluckerPolytope<P> *myPolytope = complex->getPolytope(i);
            V_TRACE(" extract stabbinglines for " << myPolytope << " at " << i);
       
   
            VisibilitySolver<P, S>::extractStabbingLines(myPolyhedron, myPolytope, mTolerance);
        }
        return VISIBLE;
    }
}
