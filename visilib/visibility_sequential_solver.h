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
            // std::count << " polytope  " << aPolytope << " is blocked " << std::endl;
            // std::count << "does not containsreallines: " << aPolytope << std::endl;
            return BLOCKED;
        }
        // std::count << "prepare to classify " << aPolytope << " against " << edges.size() <<" edges" << std::endl;
        for (auto myVisibilitySilhouetteEdge : edges)
        {            
            size_t myPolyhedronFace = myVisibilitySilhouetteEdge.mHyperPlaneIndex;
            V_ASSERT(myPolyhedronFace != 0);
            P myHyperplane = aPolyhedron->get(myPolyhedronFace);
            // std::count << "Hyperplane: " << myVisibilitySilhouetteEdge.mHyperPlaneIndex << ", " << myHyperplane << std::endl;
            auto myResult = MathPredicates::getRelativePosition(aPolytope, aPolyhedron, myHyperplane, mTolerance);
            // std::count << " classification of " << aPolytope << " against facet " << myPolyhedronFace << ":" <<myResult << std::endl;
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
         // std::count << "Preparing to split polytopes with " << edges.size() << " silhouette edges" << std::endl;
         
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
           // std::count << "Preparing to split " << workingArray.size() << " polytopes" << std::endl;
                  
             
            GeometryPositionType myResult = ON_UNKNOWN_POSITION;
            for (int i = workingArray.size() - 1; i >= 0; i--)
            {
                PluckerPolytope<P> *myPolytopeToSplit = complex->getPolytope(workingArray[i]);
                // std::count << "Split  " << myPolytopeToSplit << std::endl;
         
                PluckerPolytope<P> *myPolytopeLeft = new PluckerPolytope<P>();
                PluckerPolytope<P> *myPolytopeRight = new PluckerPolytope<P>();
                // std::count << "Creating myPolytopeLeft: " << myPolytopeLeft << std::endl;
                // std::count << "Creating myPolytopeRight: " << myPolytopeRight << std::endl;
                {
                    HelperScopedTimer timer(VisibilitySolver<P, S>::mQuery->getStatistic(), POLYTOPE_SPLIT);
                    VisibilitySolver<P, S>::mQuery->getStatistic()->inc(POLYTOPE_SPLIT_COUNT);

                    myResult = PluckerPolytopeSplitter<P, S>::split(myPolyhedron, myHyperplane, myPolytopeToSplit, myPolytopeLeft, myPolytopeRight, myPolyhedronFace, mNormalization, mTolerance);

                    if (VisibilitySolver<P, S>::mQuery->getStatistic()->get(POLYTOPE_SPLIT_COUNT) % 10000 == 0)
                    {
                        VisibilitySolver<P, S>::mQuery->getStatistic()->displayCounts();
                        // std::count << std::endl;
                    }
                }

                  // std::count << "result:  " << myResult << std::endl;
         
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
                    // std::count << "toDelete.insert:" << myPolytopeToSplit << " at index " << workingArray[i] <<  std::endl;
                    
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
                    // std::count << "Delete left and right:" << myPolytopeLeft << " " << myPolytopeRight <<  std::endl;
                  
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

          // std::count << "Preparing to classify " << complex->getPolytopeCount() << " polytopes" << std::endl;
         
        for (int polytopeIndex = 0; polytopeIndex < complex->getPolytopeCount(); polytopeIndex++)
        {
            PluckerPolytope<P> *myPolytopeToClassify = complex->getPolytope(polytopeIndex);
            // std::count << "classify " << myPolytopeToClassify << " at index " << polytopeIndex <<  std::endl;
 
           PolytopePositionType position = classifyRelativeToSilhouette(myPolytopeToClassify, myPolyhedron, edges);

            if (position == INTERSECTS)
            {
                // std::count << "adding to intersectedPolytopes "  << myPolytopeToClassify << " at " <<  polytopeIndex << std::endl;
                intersectedPolytopes.push_back(polytopeIndex);
            }
            else if (position == BLOCKED)
            {
               // std::count << "classifyRelativeToSilhouette is blocked: delete "  << myPolytopeToClassify << " at " <<  polytopeIndex << std::endl;
               

                delete myPolytopeToClassify;
                complex->setPolytope(polytopeIndex, NULL);
            }
            else 
            {
                V_ASSERT(position == OUTSIDE);
                 // std::count << "classifyRelativeToSilhouette is outside - skip: "  << myPolytopeToClassify << " at " <<  polytopeIndex << std::endl;
             
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
               // std::count << "treatIntersectedPolytopes : no intersectedpolytopes - done " << std::endl;
         
        
            return;
        }
        // std::count << "treatIntersectedPolytopes: prepare to treat " <<  intersectedPolytopes.size() << " intersected polytopes" << std::endl;
     
        for (int intersectedIndex = intersectedPolytopes.size() - 1; intersectedIndex >= 0; intersectedIndex--)
        {
            
            size_t originalPolytopeCount = complex->getPolytopeCount();
             // std::count << "originalPolytopeCount stored as  " <<  originalPolytopeCount << std::endl;
     
            std::unordered_set<PluckerPolytope<P> *> toDelete;
            std::vector<int> workingArray;
            int originalPolytopeIndex = intersectedPolytopes[intersectedIndex];

            workingArray.push_back(originalPolytopeIndex);


            // std::count << "before splitWithAllSilhouetteEdgesm working array is:  " <<  workingArray.size() << std::endl;

            splitWithAllSilhouetteEdges(workingArray, toDelete, edges);

            bool atLeastOneNewPolytopeIsBlocked = false;
     
            // std::count << "after splitWithAllSilhouetteEdges, working array is: " <<  workingArray.size() << std::endl;

            for (int i = 0; i < workingArray.size(); i++)
            {
                int myPolytopeToBlockIndex = workingArray[i];
                PluckerPolytope<P> *myPolytopeToBlock = complex->getPolytope(myPolytopeToBlockIndex);

                // std::count << "computeEdgesIntersectingQuadric for polytope " <<  myPolytopeToBlock << " at index " << myPolytopeToBlockIndex << std::endl;

                myPolytopeToBlock->computeEdgesIntersectingQuadric(myPolyhedron, mTolerance);
               
                 // std::count << "classifyRelativeToSilhouette for polytope " <<  myPolytopeToBlock << " at index " << myPolytopeToBlockIndex << std::endl;

                PolytopePositionType position = classifyRelativeToSilhouette(myPolytopeToBlock, myPolyhedron, edges);

                //V_ASSERT(position != INTERSECTS);
                if (position == BLOCKED)
                {

                 // std::count <<myPolytopeToBlock <<  " isBlocked" << std::endl;

                    if (toDelete.find(myPolytopeToBlock) != toDelete.end())
                    {
                        // std::count << "Double delete:" << myPolytopeToBlock << " " << myPolytopeToBlockIndex <<  std::endl;
                    }
                     // std::count << "toDelete.insert:" << myPolytopeToBlock << " at index " << workingArray[i] <<  std::endl;
                  
                    toDelete.insert(myPolytopeToBlock);
                    atLeastOneNewPolytopeIsBlocked = true;
                    complex->setPolytope(myPolytopeToBlockIndex, NULL);
                    workingArray[i] = -1;
                }
            }
            
            if (atLeastOneNewPolytopeIsBlocked)
            {
                 // std::count << "atLeastOneNewPolytopeIsBlocked" << std::endl;

               
                 // std::count << "set polytope at index " << originalPolytopeIndex<< " is NULL " << std::endl;

                 // std::count << "Start iterating over polytopescount to set the deleted one to NULL " << std::endl;

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
                    // std::count << "delete polytopes from delete: " << polytope<< std::endl;

                    delete polytope;
                }
            }
            else
            {
                // std::count << "!atLeastOneNewPolytopeIsBlocked: revert change " << std::endl;

                // std::count << " complex->getPolytopeCount() : " << complex->getPolytopeCount()  << " and  originalPolytopeCount:" << originalPolytopeCount << std::endl;

                while (complex->getPolytopeCount() != originalPolytopeCount)
                {
                    // std::count << "removelast" << std::endl;

                    complex->removeLast();
                     // std::count << " complex->getPolytopeCount() : " << complex->getPolytopeCount()  << " and  originalPolytopeCount:" << originalPolytopeCount << std::endl;

         
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
        // std::count << "End of procedure:" << complex->getPolytopeCount() << " left" << std::endl;
        if (complex->getPolytopeCount() == 0)
        {
             // std::count << "HIDDEN: done" << std::endl;
       
            return HIDDEN;
        }
         // std::count << "visible: prepare to extract stabbinglines"  << std::endl;
       
        for (int i = 0; i < complex->getPolytopeCount(); i++)
        {
            PluckerPolytope<P> *myPolytope = complex->getPolytope(i);
            // std::count << " extract stabbinglines for " << myPolytope << " at " << i << std::endl;
       
   
            VisibilitySolver<P, S>::extractStabbingLines(myPolyhedron, myPolytope, mTolerance);
        }
        return VISIBLE;
    }
}
