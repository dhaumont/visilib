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
    VisibilityResult VisibilitySequentialSolver<P, S>::resolve()
    {
        VisibilityResult myGlobalResult = UNKNOWN;
        Silhouette *mySilhouette = nullptr;

        bool hasEdge = false;
        bool isVisible = false;
        std::string occlusionTreeNodeSymbol;
        PluckerPolytopeComplex<P> *complex = VisibilitySolver<P, S>::mQuery->getComplex();
        PluckerPolyhedron<P> *myPolyhedron = reinterpret_cast<PluckerPolyhedron<P> *>(VisibilitySolver<P, S>::mQuery->getComplex()->getPolyhedron());

        const std::unordered_set<Silhouette *> &mySilhouettes = VisibilitySolver<P, S>::mQuery->getSilhouettes();

        for (auto iter = mySilhouettes.begin(); iter != mySilhouettes.end(); iter++)
        {
            Silhouette *mySilhouette = (*iter);
            auto &edges = mySilhouette->getEdges();
            std::vector<int> intersectedPolytopes;

            for (int polytopeIndex = 0; polytopeIndex < complex->getPolytopeCount(); polytopeIndex++)
            {
                PluckerPolytope<P> *myPolytopeToClassify = complex->getPolytope(polytopeIndex);
                bool isBlocked = true;
                bool isIntersected = false;

                for (size_t silhouetteEdgeIndex = 0; silhouetteEdgeIndex < edges.size(); silhouetteEdgeIndex++)
                {
                    GeometryPositionType myResult = ON_UNKNOWN_POSITION;
                    SilhouetteEdge &myVisibilitySilhouetteEdge = mySilhouette->getEdge(silhouetteEdgeIndex);
                    SilhouetteMeshFace *face = myVisibilitySilhouetteEdge.mFace;

                    size_t myPolyhedronFace = myVisibilitySilhouetteEdge.mHyperPlaneIndex;

                    if (myPolyhedronFace == 0)
                    {
                        MathVector2i edge = face->getEdge(myVisibilitySilhouetteEdge.mEdgeIndex);

                        MathVector3d a = convert<MathVector3d>(face->getVertex(edge.x));
                        MathVector3d b = convert<MathVector3d>(face->getVertex(edge.y));

                        P myHyperplane(a, b);
                        if (mNormalization)
                        {
                            myHyperplane = myHyperplane.getNormalized();
                        }
                        myPolyhedronFace = myPolyhedron->add(myHyperplane, ON_BOUNDARY, mNormalization, mTolerance);
                        myVisibilitySilhouetteEdge.mHyperPlaneIndex = myPolyhedronFace;
                    }

                    P myHyperplane = myPolyhedron->get(myPolyhedronFace);

                    myResult = MathPredicates::getRelativePosition(myPolytopeToClassify, myPolyhedron, myHyperplane, mTolerance);

                    if (myResult == ON_BOTH_SIDES)
                    {
                        isIntersected = true;
                    }
                    else if (myResult == ON_NEGATIVE_SIDE)
                    {
                        isBlocked = false;
                    }
                }

                if (isIntersected)
                {
                    intersectedPolytopes.push_back(polytopeIndex);
                }
                else if (isBlocked)
                {
                    delete myPolytopeToClassify;
                    complex->setPolytope(polytopeIndex, NULL);
                }
            }
            size_t originalPolytopeCount = complex->getPolytopeCount();
          
            for (int intersectedIndex = intersectedPolytopes.size() - 1; intersectedIndex >= 0; intersectedIndex--)
            {
                std::vector<PluckerPolytope<P> *> toDelete;
                std::vector<int> workingArray;

                int originalPolytopeIndex = intersectedPolytopes[intersectedIndex];
                workingArray.push_back(originalPolytopeIndex);
                for (size_t silhouetteEdgeIndex = 0; silhouetteEdgeIndex < edges.size(); silhouetteEdgeIndex++)
                {
                    SilhouetteEdge &myVisibilitySilhouetteEdge = mySilhouette->getEdge(silhouetteEdgeIndex);
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

                    GeometryPositionType myResult = ON_UNKNOWN_POSITION;
                    for (int i = workingArray.size() - 1; i >= 0; i--)
                    {
                        PluckerPolytope<P> *myPolytopeToSplit = complex->getPolytope(workingArray[i]);
                        PluckerPolytope<P> *myPolytopeLeft = new PluckerPolytope<P>();
                        PluckerPolytope<P> *myPolytopeRight = new PluckerPolytope<P>();

                        {
                            HelperScopedTimer timer(VisibilitySolver<P, S>::mQuery->getStatistic(), POLYTOPE_SPLIT);
                            VisibilitySolver<P, S>::mQuery->getStatistic()->inc(POLYTOPE_SPLIT_COUNT);

                            myResult = PluckerPolytopeSplitter<P, S>::split(myPolyhedron, myHyperplane, myPolytopeToSplit, myPolytopeLeft, myPolytopeRight, myPolyhedronFace, mNormalization, mTolerance);

                            if (VisibilitySolver<P, S>::mQuery->getStatistic()->get(POLYTOPE_SPLIT_COUNT) % 10000 == 0)
                            {
                                VisibilitySolver<P, S>::mQuery->getStatistic()->displayCounts();
                                std::cout << std::endl;
                            }
                        }

                        if (myResult == ON_BOTH_SIDES)
                        {
                            V_ASSERT(MathPredicates::getRelativePosition(myPolytopeLeft, myPolyhedron, myHyperplane, mTolerance) == ON_NEGATIVE_SIDE);
                            V_ASSERT(MathPredicates::getRelativePosition(myPolytopeRight, myPolyhedron, myHyperplane, mTolerance) == ON_POSITIVE_SIDE);
                            size_t size_before = complex->getPolytopeCount();
                            complex->appendPolytope(myPolytopeLeft);
                            complex->appendPolytope(myPolytopeRight);
                            toDelete.push_back(myPolytopeToSplit);
                            workingArray[i] = size_before;
                            workingArray.push_back(size_before + 1);
                            V_ASSERT(complex->getPolytope(workingArray[i]) == myPolytopeLeft);
                            V_ASSERT(complex->getPolytope(workingArray[workingArray.size()-1]) == myPolytopeRight);
                        }
                        else
                        {
                            V_ASSERT(MathPredicates::getRelativePosition(myPolytopeToSplit, myPolyhedron, myHyperplane, mTolerance) == myResult);
                            delete myPolytopeLeft;
                            myPolytopeLeft = nullptr;
                            delete myPolytopeRight;
                            myPolytopeRight = nullptr;
                        }
                    }
                }
                bool atLeastOnePolytopeIsBlocked = false;
                for (int i = 0; i < workingArray.size(); i++)
                {
                    int myPolytopeToBlockIndex = workingArray[i];
                    PluckerPolytope<P> *myPolytopeToBlock = complex->getPolytope(myPolytopeToBlockIndex);

                    myPolytopeToBlock->computeEdgesIntersectingQuadric(myPolyhedron, mTolerance);
                    bool isValid = myPolytopeToBlock->containsRealLines();
                    bool isBlocked = true;
                    if (isValid)
                    {
                        for (size_t silhouetteEdgeIndex = 0; silhouetteEdgeIndex < edges.size() && isBlocked; silhouetteEdgeIndex++)
                        {
                            GeometryPositionType myResult = ON_UNKNOWN_POSITION;
                            SilhouetteEdge &myVisibilitySilhouetteEdge = mySilhouette->getEdge(silhouetteEdgeIndex);
                            size_t myPolyhedronFace = myVisibilitySilhouetteEdge.mHyperPlaneIndex;
                            P myHyperplane = myPolyhedron->get(myPolyhedronFace);

                            myResult = MathPredicates::getRelativePosition(myPolytopeToBlock, myPolyhedron, myHyperplane, mTolerance);

                            if (myResult == ON_BOTH_SIDES)
                            {
                                assert(0);
                                std::cout << "ERROR: " << i << "/" << complex->getPolytopeCount() << endl;
                            }
                            if (myResult == ON_NEGATIVE_SIDE)
                            {
                                isBlocked = false;
                            }
                        }
                    }

                    if (isBlocked)
                    {
                        toDelete.push_back(myPolytopeToBlock);
                        atLeastOnePolytopeIsBlocked = true;
                        complex->setPolytope(myPolytopeToBlockIndex, NULL);
                        workingArray[i] = -1;
                    }
                }

                if (atLeastOnePolytopeIsBlocked)
                {
                    complex->setPolytope(originalPolytopeIndex, NULL);
                    for (int i = 0; i < toDelete.size(); i++)
                    {
                        delete toDelete[i];
                    }
                }
                else
                {
                    while (complex->getPolytopeCount() != originalPolytopeCount)
                    {
                        complex->removeLast();
                    }                    
                }
            }

            complex->repack();
        }
        if (complex->getPolytopeCount() == 0)
        {
            return HIDDEN;
        }
        for (int i = complex->getPolytopeCount() - 1; i >= 0; i--)
        {
            PluckerPolytope<P> *myPolytope = complex->getPolytope(i);
            if (myPolytope->containsRealLines())
            {
                VisibilitySolver<P, S>::extractStabbingLines(myPolyhedron, myPolytope, mTolerance);
            }
        }
        return VISIBLE;
    }
}
