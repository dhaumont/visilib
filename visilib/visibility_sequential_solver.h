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
#include "geometry_ray.h"
#include "silhouette_mesh_face.h"
#include "math_geometry.h"
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
#include <set>
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
        /* Option B blocked predicate: cast the polytope's representative line
        (converted back to 3D) against the faces of the given silhouette only.
        Returns true when the line is occluded by that patch, false otherwise.
        A polytope with no real lines is always occluded. */

        bool isBlockedBySilhouette(PluckerPolytope<P> *aPolytope,
                                   PluckerPolyhedron<P> *aPolyhedron,
                                   const Silhouette *aSilhouette);

        PolytopePositionType classifyRelativeToSilhouette(PluckerPolytope<P> *myPolytopeToClassify,
                                          PluckerPolyhedron<P> *myPolyhedron,
                                          const std::vector<SilhouetteEdge> &edges,
                                          const Silhouette *aSilhouette
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
            std::vector<int> &intersectedPolytopes,
            const Silhouette *aSilhouette);

        void treatIntersectedPolytopes(
            const std::vector<SilhouetteEdge> &edges,
            PluckerPolytopeComplex<P> *complex,
            PluckerPolyhedron<P> *myPolyhedron,
            const std::vector<int> &intersectedPolytopes,
            const Silhouette *aSilhouette);

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
            if (VisibilitySolver<P, S>::isTimedOut())
            {
                return;
            }

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
    bool VisibilitySequentialSolver<P, S>::isBlockedBySilhouette(PluckerPolytope<P> *aPolytope,
                                                                  PluckerPolyhedron<P> *aPolyhedron,
                                                                  const Silhouette *aSilhouette)
    {
        if (!aPolytope->containsRealLines())
        {
            return true;
        }

        P repLine = MathGeometry::computeRepresentativeLine<P>(aPolytope, aPolyhedron, mTolerance);

        const MathPlane3d &plane0 = VisibilitySolver<P, S>::mQuery->getQueryPolygon(0)->getPlane();
        const MathPlane3d &plane1 = VisibilitySolver<P, S>::mQuery->getQueryPolygon(1)->getPlane();

        std::pair<MathVector3d, MathVector3d> seg = MathGeometry::getBackTo3D(repLine, plane0, plane1);

        MathVector3d dir = seg.second - seg.first;
        const double len = dir.getNorm();
        if (len <= 0.0)
        {
            return false;
        }
        dir *= (1.0 / len);

        MathVector3f start((float)seg.first.x, (float)seg.first.y, (float)seg.first.z);
        MathVector3f direction((float)dir.x, (float)dir.y, (float)dir.z);
        GeometryRay ray(start, direction);

        /* Test the representative line against all faces of this silhouette's
        occluder geometry, not just the faces in the silhouette patch.  The patch
        only contains faces that project into the source-to-target shaft; a rep
        line might exit through an adjacent face that is not in the patch (e.g. the
        far side of a closed convex occluder).

        getMeshFaces() returns the full face array for this occluder (the same
        buffer used by the query's silhouette processor), which is the correct
        geometric reference independent of which occluder set accessor we use. */

        const std::vector<SilhouetteMeshFace> &allFaces = aSilhouette->getMeshFaces();

        for (const SilhouetteMeshFace &face : allFaces)
        {
            if (MathGeometry::hitsTriangle<float>(ray, face.getVertex(0), face.getVertex(1), face.getVertex(2)))
            {
                return true;
            }
        }
        return false;
    }

    template <class P, class S>
    PolytopePositionType VisibilitySequentialSolver<P, S>::classifyRelativeToSilhouette(PluckerPolytope<P> *aPolytope,
                                                                        PluckerPolyhedron<P> *aPolyhedron,
                                                                        const std::vector<SilhouetteEdge> &edges,
                                                                        const Silhouette *aSilhouette)
    {
        if (!aPolytope->containsRealLines())
        {
            V_TRACE("polytope " << aPolytope << " has no real lines: BLOCKED");
            return BLOCKED;
        }

        /* First pass: check whether the polytope straddles any silhouette hyperplane.
        A cell that crosses a hyperplane must be split before we can classify it —
        only a fully one-sided cell has a well-defined representative line. */

        bool straddled = false;
        V_TRACE("classify " << aPolytope << " against " << edges.size() << " edges");
        for (auto &myEdge : edges)
        {
            size_t myPolyhedronFace = myEdge.mHyperPlaneIndex;
            V_ASSERT(myPolyhedronFace != 0);
            P myHyperplane = aPolyhedron->get(myPolyhedronFace);
            auto myResult = MathPredicates::getRelativePosition(aPolytope, aPolyhedron, myHyperplane, mTolerance);
            V_TRACE("facet " << myPolyhedronFace << ": " << myResult);
            if (myResult == ON_BOTH_SIDES)
            {
                straddled = true;
                break;
            }
        }

        if (straddled)
        {
            V_TRACE("polytope " << aPolytope << " straddles a hyperplane: INTERSECTS");
            return INTERSECTS;
        }

        /* Second pass (Option B): the polytope is on one side of every hyperplane.
        Cast its representative line against the silhouette faces to decide. */

        if (isBlockedBySilhouette(aPolytope, aPolyhedron, aSilhouette))
        {
            V_TRACE("polytope " << aPolytope << " representative line hits occluder: BLOCKED");
            return BLOCKED;
        }

        V_TRACE("polytope " << aPolytope << " clears occluder: OUTSIDE");
        return OUTSIDE;
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
            /* Every loop of this solver gives up on an over budget query, and only resolve()
            draws conclusions, which is where the budget is checked again. */
            if (VisibilitySolver<P, S>::isTimedOut())
            {
                return;
            }

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

                    if (VisibilitySolver<P, S>::mQuery->getStatistic()->hasToDisplaySplitCount())
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
                        /* F6: the split produced two empty children — the parent is
                        fully consumed by the hyperplane.  Remove its slot from the
                        working array so the caller never touches it again. */
                        V_TRACE("both children empty after split: parent consumed, removing from workingArray");
                        workingArray.erase(workingArray.begin() + i);
                        i--;
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
        std::vector<int> &intersectedPolytopes,
        const Silhouette *aSilhouette)
    {
        V_TRACE("Preparing to classify " << complex->getPolytopeCount() << " polytopes");

        for (int polytopeIndex = 0; polytopeIndex < complex->getPolytopeCount(); polytopeIndex++)
        {
            if (VisibilitySolver<P, S>::isTimedOut())
            {
                return;
            }

            PluckerPolytope<P> *myPolytopeToClassify = complex->getPolytope(polytopeIndex);
            V_TRACE("classify " << myPolytopeToClassify << " at index " << polytopeIndex);

            PolytopePositionType position = classifyRelativeToSilhouette(myPolytopeToClassify, myPolyhedron, edges, aSilhouette);

            if (position == INTERSECTS)
            {
                V_TRACE("adding to intersectedPolytopes " << myPolytopeToClassify << " at " << polytopeIndex);
                intersectedPolytopes.push_back(polytopeIndex);
            }
            else if (position == BLOCKED)
            {
                V_TRACE("classifyRelativeToSilhouette is blocked: delete " << myPolytopeToClassify << " at " << polytopeIndex);
                delete myPolytopeToClassify;
                complex->setPolytope(polytopeIndex, NULL);
            }
            else
            {
                V_ASSERT(position == OUTSIDE);
                V_TRACE("classifyRelativeToSilhouette is outside - skip: " << myPolytopeToClassify << " at " << polytopeIndex);
            }
        }
    }

    template <class P, class S>
    void VisibilitySequentialSolver<P, S>::treatIntersectedPolytopes(
        const std::vector<SilhouetteEdge> &edges,
        PluckerPolytopeComplex<P> *complex,
        PluckerPolyhedron<P> *myPolyhedron,
        const std::vector<int> &intersectedPolytopes,
        const Silhouette *aSilhouette)
    {
        if (intersectedPolytopes.empty())
        {
            V_TRACE("treatIntersectedPolytopes: no intersected polytopes - done");
            return;
        }
        V_TRACE("treatIntersectedPolytopes: prepare to treat " << intersectedPolytopes.size() << " intersected polytopes");

        for (int intersectedIndex = intersectedPolytopes.size() - 1; intersectedIndex >= 0; intersectedIndex--)
        {
            size_t originalPolytopeCount = complex->getPolytopeCount();
            V_TRACE("originalPolytopeCount stored as " << originalPolytopeCount);

            std::unordered_set<PluckerPolytope<P> *> toDelete;
            std::vector<int> workingArray;
            int originalPolytopeIndex = intersectedPolytopes[intersectedIndex];
            workingArray.push_back(originalPolytopeIndex);

            V_TRACE("before splitWithAllSilhouetteEdges, working array size: " << workingArray.size());

            splitWithAllSilhouetteEdges(workingArray, toDelete, edges);

            bool atLeastOneNewPolytopeIsBlocked = false;

            V_TRACE("after splitWithAllSilhouetteEdges, working array size: " << workingArray.size());

            for (int i = 0; i < (int)workingArray.size(); i++)
            {
                if (workingArray[i] < 0)
                {
                    continue;
                }
                int myPolytopeToBlockIndex = workingArray[i];
                PluckerPolytope<P> *myPolytopeToBlock = complex->getPolytope(myPolytopeToBlockIndex);
                if (myPolytopeToBlock == nullptr)
                {
                    continue;
                }

                V_TRACE("computeEdgesIntersectingQuadric for polytope " << myPolytopeToBlock << " at index " << myPolytopeToBlockIndex);
                myPolytopeToBlock->computeEdgesIntersectingQuadric(myPolyhedron, mTolerance);

                V_TRACE("classifyRelativeToSilhouette for polytope " << myPolytopeToBlock << " at index " << myPolytopeToBlockIndex);
                PolytopePositionType position = classifyRelativeToSilhouette(myPolytopeToBlock, myPolyhedron, edges, aSilhouette);

                if (position == BLOCKED)
                {
                    V_TRACE(myPolytopeToBlock << " isBlocked");
                    V_TRACE("toDelete.insert: " << myPolytopeToBlock << " at index " << workingArray[i]);
                    toDelete.insert(myPolytopeToBlock);
                    atLeastOneNewPolytopeIsBlocked = true;
                    complex->setPolytope(myPolytopeToBlockIndex, NULL);
                    workingArray[i] = -1;
                }
            }

            if (atLeastOneNewPolytopeIsBlocked)
            {
                V_TRACE("atLeastOneNewPolytopeIsBlocked");

                for (int i = 0; i < (int)complex->getPolytopeCount(); i++)
                {
                    PluckerPolytope<P> *myPolytope = complex->getPolytope(i, false);
                    if (toDelete.find(myPolytope) != toDelete.end())
                    {
                        complex->setPolytope(i, NULL);
                    }
                }

                for (auto polytope : toDelete)
                {
                    V_TRACE("delete polytope: " << polytope);
                    delete polytope;
                }
            }
            else
            {
                V_TRACE("!atLeastOneNewPolytopeIsBlocked: revert change");
                V_TRACE("complex->getPolytopeCount(): " << complex->getPolytopeCount() << " originalPolytopeCount: " << originalPolytopeCount);

                while (complex->getPolytopeCount() != originalPolytopeCount)
                {
                    V_TRACE("removelast");
                    complex->removeLast();
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

        std::map<size_t, std::vector<Silhouette*>> myOccluders;

        for (auto iter = mySilhouettes.begin(); iter != mySilhouettes.end(); iter++)
        {
            myOccluders[(*iter)->getGeometryId()].push_back(*iter);
        }

        for (auto occluder = myOccluders.begin(); occluder != myOccluders.end(); occluder++)
        {
            if (VisibilitySolver<P, S>::isTimedOut())
            {
                V_TRACE("FAILURE: time budget over");
                return FAILURE;
            }

            std::vector<Silhouette*> &myPatches = occluder->second;
            Silhouette *mySilhouette = myPatches[0];

            std::vector<SilhouetteEdge> myEdges;
            std::set<std::pair<SilhouetteMeshFace*, size_t>> mySeen;

            for (size_t patch = 0; patch < myPatches.size(); patch++)
            {
                for (SilhouetteEdge &myEdge : myPatches[patch]->getEdges())
                {
                    if (mySeen.insert(std::make_pair(myEdge.mFace, myEdge.mEdgeIndex)).second)
                    {
                        myEdges.push_back(myEdge);
                    }
                }
            }

            std::vector<SilhouetteEdge> &edges = myEdges;

            if (edges.empty())
            {
                V_TRACE("silhouette has no edges: apply ray-test-only classification");
                for (int polytopeIndex = 0; polytopeIndex < complex->getPolytopeCount(); polytopeIndex++)
                {
                    PluckerPolytope<P> *myPolytope = complex->getPolytope(polytopeIndex);
                    if (myPolytope != nullptr && isBlockedBySilhouette(myPolytope, myPolyhedron, mySilhouette))
                    {
                        V_TRACE("polytope " << myPolytope << " blocked by empty-edge silhouette: delete");
                        delete myPolytope;
                        complex->setPolytope(polytopeIndex, nullptr);
                    }
                }
                complex->repack();
                continue;
            }

            createMissingHyperplanes(myPolyhedron, edges);

            std::vector<int> intersectedPolytopes;
            classifyAllPolytopes(edges, complex, myPolyhedron, intersectedPolytopes, mySilhouette);

            treatIntersectedPolytopes(edges, complex, myPolyhedron, intersectedPolytopes, mySilhouette);

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
