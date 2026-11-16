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

#include <algorithm>
#include <unordered_set>

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
#include "shadow_discontinuity_tracer.h"
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <vector>

namespace visilib
{
    /** @brief One silhouette triangle, with its three edge hyperplanes calibrated.

    A plane is stored as the polyhedron face of its edge, given the sign of the side a line
    through the face centroid stands on: that line crosses the triangle, so the three signed
    planes describe the same thing whatever orientation the edges carry, and no orientation
    question is left to the tests made on them. */

    struct FaceShadow
    {
        const SilhouetteMeshFace *mFace = nullptr;
        int mPlanes[3] = { 0, 0, 0 };
    };

    /** @brief How a cell of lines stands with respect to one calibrated face.

    The crossing of the plane of a triangle carries its three edge forms alike: a line hitting
    it is strictly inside the three planes, or, taken from the back, strictly outside the three.
    A cell inside one plane and outside another therefore holds nothing that reaches the face,
    and a cell uniform on the three holds nothing but lines that do. Anything less strict, a
    cell touching a plane or resting on one vertex, settles nothing and is left unknown. */

    enum CellFacePosition
    {
        FACE_UNKNOWN = 0,   /**< @brief the cell is not strict enough on the planes for a verdict*/
        FACE_BLOCKS,        /**< @brief every line of the cell hits the face*/
        FACE_AVOIDS         /**< @brief no line of the cell hits the face*/
    };

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
        /* Option B blocked predicate: cast the representative line of the polytope
        (converted back to 3D) against the faces of the given silhouette only.
        Returns BLOCKED when the line is occluded by that patch, OUTSIDE when it clears it,
        and UNDEFINED when no line inside the polytope could be found to cast, which is decided
        only by the optimised procedure: the gravity center one is used as it stands.
        A polytope with no real lines is always occluded.
        The result is an enum where BLOCKED is nought: compare it, never test it as a condition. */

        PolytopePositionType isBlockedBySilhouette(PluckerPolytope<P> *aPolytope,
                                   PluckerPolyhedron<P> *aPolyhedron,
                                   const Silhouette *aSilhouette);

        PolytopePositionType classifyRelativeToSilhouette(PluckerPolytope<P> *myPolytopeToClassify,
                                          PluckerPolyhedron<P> *myPolyhedron,
                                          const std::vector<SilhouetteEdge> &edges,
                                          const Silhouette *aSilhouette,
                                          bool skipStraddleScan = false
                                          );

        void splitWithAllSilhouetteEdges(std::vector<int> &workingArray,
                                         std::unordered_set<PluckerPolytope<P> *> &toDelete,
                                         const std::vector<SilhouetteEdge> &edges,
                                         const std::map<const SilhouetteMeshFace *, FaceShadow> &shadows,
                                         bool &earlyBlocked);

        void createMissingHyperplanes(PluckerPolyhedron<P> *aPolyhedron,
                                      std::vector<SilhouetteEdge> &edges);

        /* The calibrated faces of a silhouette: one entry per triangle whose three edges are
        all of the pass, once their planes exist. Dense occluders only, whose arrangement is
        the one that never ends without these tests; a sparse silhouette gets an empty map and
        every caller falls back on the pass as it was. */

        void collectFaceShadows(PluckerPolyhedron<P> *aPolyhedron,
                                const std::vector<SilhouetteEdge> &edges,
                                std::map<const SilhouetteMeshFace *, FaceShadow> &shadows) const;

        bool isPolytopeStrictlyInside(PluckerPolytope<P> *aPolytope,
                                      PluckerPolyhedron<P> *aPolyhedron,
                                      int calibratedPlaneId) const;

        CellFacePosition cellAgainstFace(PluckerPolytope<P> *aPolytope,
                                         PluckerPolyhedron<P> *aPolyhedron,
                                         const FaceShadow &shadow) const;

        /* The faces of the occluder on duty that no line of the beam can reach: out of the
        edges of the pass, and out of the ray cast. A face is dropped with its three planes
        and kept with them, never one without the other: the classification concludes that a
        cell is clear of the faces it still casts, and that answer needs the planes of every
        face it keeps. A face missing from the calibrated set is kept: it was not proved of
        anything. */

        void filterUnreachableFaces(const std::map<const SilhouetteMeshFace *, FaceShadow> &shadows,
                                    PluckerPolyhedron<P> *aPolyhedron,
                                    PluckerPolytope<P> *beam,
                                    std::vector<SilhouetteEdge> &edges);

        void classifyAllPolytopes(
            const std::vector<SilhouetteEdge> &edges,
            PluckerPolytopeComplex<P> *complex,
            PluckerPolyhedron<P> *myPolyhedron,
            std::vector<int> &intersectedPolytopes,
            const Silhouette *aSilhouette);

        void treatIntersectedPolytopes(
            const std::vector<SilhouetteEdge> &edges,
            const std::map<const SilhouetteMeshFace *, FaceShadow> &shadows,
            PluckerPolytopeComplex<P> *complex,
            PluckerPolyhedron<P> *myPolyhedron,
            const std::vector<int> &intersectedPolytopes,
            const Silhouette *aSilhouette);

        bool mNormalization;
        bool mDetectApertureOnly;
        bool mOptimisedRepresentativeLine;   /**< @brief Whether the representative line of a cell is looked for inside it, or taken from the gravity center*/
        S mTolerance;

        /* The faces of the occluder on duty that no line of the beam reaches, proved against
        the beam itself before any cell was cut. They are out of the ray cast of the blocked
        predicate: a face no line of the beam touches cannot decide anything for any cell, and
        the arrangement keeps the planes of every face it does cast. */
        std::unordered_set<const SilhouetteMeshFace *> mUnreachableFaces;
    };

    template <class P, class S>
    VisibilitySequentialSolver<P, S>::VisibilitySequentialSolver(VisibilityExactQuery_<P, S> *mQuery, bool normalization, S tolerance, bool detectApertureOnly)
        : VisibilitySolver<P, S>(mQuery),
          mNormalization(normalization),
          mTolerance(tolerance),
          mDetectApertureOnly(detectApertureOnly)
    {
        /* How the representative line of a cell is found belongs to the query, not to the caller: the interface
        switches it to compare the two procedures, everything else runs the one the configuration asks for. */
        mOptimisedRepresentativeLine = mQuery->getConfiguration().optimisedRepresentativeLine;
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

                P myHyperplane = MathGeometry::createEdgeFromPoints<P>(a, b);
                if (mNormalization)
                {
                    myHyperplane = myHyperplane.getNormalized();
                }
                myPolyhedronFace = aPolyhedron->add(myHyperplane, ON_BOUNDARY, mNormalization, mTolerance);
                myVisibilitySilhouetteEdge.mHyperPlaneIndex = myPolyhedronFace;

                HyperplaneEventInfo info;
                info.mKind = HyperplaneEventInfo::SilhouetteEdge;
                info.mEdgeIndex = myVisibilitySilhouetteEdge.mEdgeIndex;
                info.mFaceIndex = face->getFaceIndex();
                info.mA = a;
                info.mB = b;
                VisibilitySolver<P, S>::mQuery->getEventRegistry().registerHyperplane(myPolyhedronFace, info);
            }
        }
    }

    template <class P, class S>
    void VisibilitySequentialSolver<P, S>::collectFaceShadows(PluckerPolyhedron<P> *aPolyhedron,
                                                              const std::vector<SilhouetteEdge> &edges,
                                                              std::map<const SilhouetteMeshFace *, FaceShadow> &shadows) const
    {
        shadows.clear();

        if (edges.size() <= 300)
        {
            return;
        }

        std::map<const SilhouetteMeshFace *, std::vector<const SilhouetteEdge *>> grouped;
        for (const SilhouetteEdge &myEdge : edges)
        {
            if (myEdge.mFace != nullptr)
            {
                grouped[myEdge.mFace].push_back(&myEdge);
            }
        }

        for (const auto &group : grouped)
        {
            /* Only a triangle whose three edges are of the pass can be calibrated: a face
            missing one of them has a plane nobody created, and the cone it draws is
            incomplete. */
            if (group.second.size() != 3)
            {
                continue;
            }

            const SilhouetteMeshFace *face = group.first;
            MathVector3f v0 = face->getVertex(0);
            MathVector3f v1 = face->getVertex(1);
            MathVector3f v2 = face->getVertex(2);
            MathVector3f normal = MathVector3f::cross(v1 - v0, v2 - v0);
            if (normal.getSquaredNorm() <= 0.0f)
            {
                continue;
            }
            MathVector3f centroid((v0.x + v1.x + v2.x) / 3.0f, (v0.y + v1.y + v2.y) / 3.0f, (v0.z + v1.z + v2.z) / 3.0f);
            MathVector3d testPoint = convert<MathVector3d>(MathVector3f(centroid.x + normal.x, centroid.y + normal.y, centroid.z + normal.z));
            MathVector3d testOrigin = convert<MathVector3d>(MathVector3f(centroid.x - normal.x, centroid.y - normal.y, centroid.z - normal.z));
            P testLine = MathGeometry::createEdgeFromPoints<P>(testOrigin, testPoint);
            if (mNormalization)
            {
                testLine = testLine.getNormalized();
            }

            FaceShadow shadow;
            shadow.mFace = face;
            bool calibrated = true;
            for (int k = 0; k < 3; k++)
            {
                const SilhouetteEdge &se = *group.second[k];
                if (se.mHyperPlaneIndex == 0)
                {
                    calibrated = false;
                    break;
                }
                const P &plane = aPolyhedron->get(se.mHyperPlaneIndex);
                S side = plane.dot(testLine);
                /* A centroid line standing exactly on a plane of its own face draws no
                inside at all: the face is left aside. The exact scalars have no equality,
                and the two strict signs say everything there is to know. */
                if (!(side > S(0)) && !(side < S(0)))
                {
                    calibrated = false;
                    break;
                }
                shadow.mPlanes[k] = side > S(0) ? (int)se.mHyperPlaneIndex : -(int)se.mHyperPlaneIndex;
            }

            if (calibrated)
            {
                shadows[face] = shadow;
            }
        }
    }

    template <class P, class S>
    bool VisibilitySequentialSolver<P, S>::isPolytopeStrictlyInside(PluckerPolytope<P> *aPolytope,
                                                                    PluckerPolyhedron<P> *aPolyhedron,
                                                                    int calibratedPlaneId) const
    {
        if (aPolytope == nullptr)
        {
            return false;
        }

        const std::unordered_set<size_t> &vertices = aPolytope->getVertices();
        if (vertices.empty())
        {
            return false;
        }

        const size_t plane = calibratedPlaneId < 0 ? (size_t)-calibratedPlaneId : (size_t)calibratedPlaneId;
        const P &hyperplane = aPolyhedron->get(plane);
        const GeometryPositionType wanted = calibratedPlaneId > 0 ? ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE;

        for (size_t v : vertices)
        {
            if (MathPredicates::getRelativePosition(hyperplane.dot(aPolyhedron->get(v)), mTolerance) != wanted)
            {
                return false;
            }
        }
        return true;
    }

    template <class P, class S>
    CellFacePosition VisibilitySequentialSolver<P, S>::cellAgainstFace(PluckerPolytope<P> *aPolytope,
                                                                       PluckerPolyhedron<P> *aPolyhedron,
                                                                       const FaceShadow &shadow) const
    {
        int side[3] = { 0, 0, 0 };
        for (int k = 0; k < 3; k++)
        {
            if (isPolytopeStrictlyInside(aPolytope, aPolyhedron, shadow.mPlanes[k]))
            {
                side[k] = 1;
            }
            else if (isPolytopeStrictlyInside(aPolytope, aPolyhedron, -shadow.mPlanes[k]))
            {
                side[k] = -1;
            }
        }

        for (int k = 0; k < 3; k++)
        {
            for (int j = k + 1; j < 3; j++)
            {
                if (side[k] != 0 && side[j] != 0 && side[k] != side[j])
                {
                    return FACE_AVOIDS;
                }
            }
        }

        if (side[0] != 0 && side[1] != 0 && side[2] != 0)
        {
            return FACE_BLOCKS;
        }
        return FACE_UNKNOWN;
    }

    template <class P, class S>
    void VisibilitySequentialSolver<P, S>::filterUnreachableFaces(const std::map<const SilhouetteMeshFace *, FaceShadow> &shadows,
                                                                  PluckerPolyhedron<P> *aPolyhedron,
                                                                  PluckerPolytope<P> *beam,
                                                                  std::vector<SilhouetteEdge> &edges)
    {
        mUnreachableFaces.clear();

        if (shadows.empty() || beam == nullptr)
        {
            return;
        }

        for (const auto &shadow : shadows)
        {
            if (cellAgainstFace(beam, aPolyhedron, shadow.second) == FACE_AVOIDS)
            {
                mUnreachableFaces.insert(shadow.first);
            }
        }

        if (mUnreachableFaces.empty())
        {
            return;
        }

        const size_t edgesBefore = edges.size();
        edges.erase(std::remove_if(edges.begin(), edges.end(),
                                   [this](const SilhouetteEdge &anEdge)
                                   {
                                       return mUnreachableFaces.find(anEdge.mFace) != mUnreachableFaces.end();
                                   }),
                    edges.end());

        HelperStatisticCollector *statistic = VisibilitySolver<P, S>::mQuery->getStatistic();
        statistic->add(SILHOUETTE_FACE_REJECTED_COUNT, mUnreachableFaces.size());
        statistic->add(SILHOUETTE_EDGE_DROPPED_COUNT, edgesBefore - edges.size());

        V_TRACE("beam reaches none of " << mUnreachableFaces.size() << " silhouette faces: "
                                        << edgesBefore << " edges down to " << edges.size());
    }

    template <class P, class S>
    PolytopePositionType VisibilitySequentialSolver<P, S>::isBlockedBySilhouette(PluckerPolytope<P> *aPolytope,
                                                                  PluckerPolyhedron<P> *aPolyhedron,
                                                                  const Silhouette *aSilhouette)
    {
        if (!aPolytope->containsRealLines())
        {
            return BLOCKED;
        }

        const RepresentativeLine<P> representative = MathGeometry::representativeLine<P, S>(aPolytope, aPolyhedron, mTolerance, mOptimisedRepresentativeLine,
                                                            VisibilitySolver<P, S>::mQuery->getStatistic());

        /* A line that lies outside the cell tells nothing about it: it may graze an occluder the cell itself does
        not see, or miss one it does. Asked for a certificate, the search says it does not know and the cell is
        split further instead of concluded from. */
        if (mOptimisedRepresentativeLine && !representative.insideCell && !representative.boundary)
        {
            V_TRACE("polytope " << aPolytope << " has no verified representative line (margin " << representative.margin << ")");
            return UNDEFINED;
        }

        P repLine = representative.line;

        const MathPlane3d &plane0 = VisibilitySolver<P, S>::mQuery->getQueryPolygon(0)->getPlane();
        const MathPlane3d &plane1 = VisibilitySolver<P, S>::mQuery->getQueryPolygon(1)->getPlane();

        std::pair<MathVector3d, MathVector3d> seg = MathGeometry::getBackTo3D(repLine, plane0, plane1);

        MathVector3d dir = seg.second - seg.first;
        const double len = dir.getNorm();
        if (len <= 0.0)
        {
            return mOptimisedRepresentativeLine ? UNDEFINED : OUTSIDE;
        }
        dir *= (1.0 / len);

        MathVector3f start((float)seg.first.x, (float)seg.first.y, (float)seg.first.z);
        MathVector3f direction((float)dir.x, (float)dir.y, (float)dir.z);
        GeometryRay ray(start, direction);

        const std::vector<SilhouetteMeshFace> &allFaces = aSilhouette->getMeshFaces();

        for (const SilhouetteMeshFace &face : allFaces)
        {
            /* A face the beam cannot reach is out of the question asked here, and its own
            edges are still in the arrangement only as long as another face needs them. */
            if (!mUnreachableFaces.empty() && mUnreachableFaces.find(&face) != mUnreachableFaces.end())
            {
                continue;
            }
            if (MathGeometry::hitsTriangle<float>(ray, face.getVertex(0), face.getVertex(1), face.getVertex(2)))
            {
                return BLOCKED;
            }
        }
        return OUTSIDE;
    }

    template <class P, class S>
    PolytopePositionType VisibilitySequentialSolver<P, S>::classifyRelativeToSilhouette(PluckerPolytope<P> *aPolytope,
                                                                        PluckerPolyhedron<P> *aPolyhedron,
                                                                        const std::vector<SilhouetteEdge> &edges,
                                                                        const Silhouette *aSilhouette,
                                                                        bool skipStraddleScan)
    {
        if (!aPolytope->containsRealLines())
        {
            V_TRACE("polytope " << aPolytope << " has no real lines: BLOCKED");
            return BLOCKED;
        }

        /* First pass: check whether the polytope straddles any silhouette hyperplane.
        A cell that crosses a hyperplane must be split before we can classify it —
        only a fully one-sided cell has a well-defined representative line.

        The split pass leaves no cell straddling the occluder's planes: a child is born on
        one side of the plane that cut it, inherits from its parent the sides it was already
        on, and the pass carries every survivor through the planes still to come. Asked from
        just out of a pass run to its end, the scan below re-reads what the splitting proved
        plane by plane, and is skipped. */

        bool straddled = false;
        V_TRACE("classify " << aPolytope << " against " << edges.size() << " edges");
        if (!skipStraddleScan)
        {
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
        }

        if (straddled)
        {
            V_TRACE("polytope " << aPolytope << " straddles a hyperplane: INTERSECTS");
            return INTERSECTS;
        }

        /* Second pass (Option B): the polytope is on one side of every hyperplane.
        Cast its representative line against the silhouette faces to decide. */

        const PolytopePositionType position = isBlockedBySilhouette(aPolytope, aPolyhedron, aSilhouette);
        V_TRACE("polytope " << aPolytope << " representative line decides " << position);
        return position;
    }

    template <class P, class S>
    void VisibilitySequentialSolver<P, S>::splitWithAllSilhouetteEdges(std::vector<int> &workingArray,
                                                                       std::unordered_set<PluckerPolytope<P> *> &toDelete,
                                                                       const std::vector<SilhouetteEdge> &edges,
                                                                       const std::map<const SilhouetteMeshFace *, FaceShadow> &shadows,
                                                                       bool &earlyBlocked)
    {
        PluckerPolytopeComplex<P> *complex = VisibilitySolver<P, S>::mQuery->getComplex();
        PluckerPolyhedron<P> *myPolyhedron = reinterpret_cast<PluckerPolyhedron<P> *>(VisibilitySolver<P, S>::mQuery->getComplex()->getPolyhedron());
         V_TRACE("Preparing to split polytopes with " << edges.size() << " silhouette edges");

        earlyBlocked = false;

        /* The pass walks the edges grouped by face when the shadows are there: the three
        planes of a face then come one after the other, and a cell dies at the third of them
        instead of living through the pass to be buried by the classification one by one.
        Without shadows, the edges are taken as they were given. */

        std::vector<size_t> orderedEdges(edges.size());
        for (size_t e = 0; e < orderedEdges.size(); e++)
        {
            orderedEdges[e] = e;
        }

        if (!shadows.empty())
        {
            std::sort(orderedEdges.begin(), orderedEdges.end(),
                      [&edges](size_t a, size_t b)
                      {
                          if (edges[a].mFace != edges[b].mFace)
                          {
                              return std::less<const SilhouetteMeshFace *>()(edges[a].mFace, edges[b].mFace);
                          }
                          return a < b;
                      });
        }

        for (size_t orderedIndex = 0; orderedIndex < orderedEdges.size(); orderedIndex++)
        {
             if (VisibilitySolver<P, S>::isTimedOut())
            {
                return;
            }
            const SilhouetteEdge &myVisibilitySilhouetteEdge = edges[orderedEdges[orderedIndex]];
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
                if (VisibilitySolver<P, S>::isTimedOut())
                {
                    return;
                }

                if (workingArray[i] < 0)
                {
                    /* A slot retired by the pruning below: the polytope it held is gone and
                    the remaining edges have nothing left to cut there. */
                    continue;
                }

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

                    /* The cell was not split by this edge: it lies on one side of it, and one
                    side of a plane of its face is one side of the three asked for. Uniform on
                    them, every line it holds hits that face: it is blocked, and it goes now,
                    before the edges still to come have carved anything else out of it. It is
                    the verdict the classification after the pass would reach, taken as soon as
                    it is provable — with a dense occluder the difference between an answer and
                    an arrangement that never ends. */

                    if (myResult != ON_BOUNDARY && !shadows.empty())
                    {
                        auto shadow = shadows.find(face);
                        if (shadow != shadows.end() && cellAgainstFace(myPolytopeToSplit, myPolyhedron, shadow->second) == FACE_BLOCKS)
                        {
                            V_TRACE("polytope at slot " << workingArray[i] << " lies in the shadow of one face: blocked");
                            VisibilitySolver<P, S>::mQuery->getStatistic()->inc(POLYTOPE_SHADOW_BLOCKED_COUNT);
                            toDelete.insert(myPolytopeToSplit);
                            complex->setPolytope(workingArray[i], NULL);
                            workingArray[i] = -1;
                            earlyBlocked = true;
                        }
                    }
                  
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

            /* A cell that could not be decided is treated as one the silhouette cuts across: it gets split, which
            costs a split and loses nothing. Concluding from a line that is not inside the cell is what loses. */
            if (position == INTERSECTS || position == UNDEFINED)
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
        const std::map<const SilhouetteMeshFace *, FaceShadow> &shadows,
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
           if (VisibilitySolver<P, S>::isTimedOut())
            {
                return;
            }
            size_t originalPolytopeCount = complex->getPolytopeCount();
            V_TRACE("originalPolytopeCount stored as " << originalPolytopeCount);

            std::unordered_set<PluckerPolytope<P> *> toDelete;
            std::vector<int> workingArray;
            int originalPolytopeIndex = intersectedPolytopes[intersectedIndex];
            workingArray.push_back(originalPolytopeIndex);

            V_TRACE("before splitWithAllSilhouetteEdges, working array size: " << workingArray.size());

            bool earlyBlocked = false;
            splitWithAllSilhouetteEdges(workingArray, toDelete, edges, shadows, earlyBlocked);

            /* The survivors of a pass run to its end are on one side of every plane of the
            occluder: the classification below asks them that question no more. A pass given
            up on mid-way has not visited the planes after the one it stopped at, and its
            cells are scanned as any others would be. */
            const bool splitPassRanToTheEnd = !VisibilitySolver<P, S>::isTimedOut();

            /* A cell the pass buried in a face shadow is a new blocked polytope of this
            occluder: the verdict is taken, and the pass counts as one that found blockage. */
            bool atLeastOneNewPolytopeIsBlocked = earlyBlocked;

            V_TRACE("after splitWithAllSilhouetteEdges, working array size: " << workingArray.size());

            for (int i = 0; i < (int)workingArray.size(); i++)
            {
                if (VisibilitySolver<P, S>::isTimedOut())
                {
                    return;
                }

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
                PolytopePositionType position = classifyRelativeToSilhouette(myPolytopeToBlock, myPolyhedron, edges, aSilhouette, splitPassRanToTheEnd);

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

        const std::vector<Silhouette *> &mySilhouettes = VisibilitySolver<P, S>::mQuery->getSilhouettes();

        std::map<size_t, std::vector<Silhouette*>> myOccluders;

        for (auto iter = mySilhouettes.begin(); iter != mySilhouettes.end(); iter++)
        {
            myOccluders[(*iter)->getGeometryId()].push_back(*iter);
        }

        /* The beam, kept whole: the cell the complex starts from is the set of lines through
        the two query polygons, and every cell the pass will cut is a piece of it. Held aside,
        it answers what a silhouette face can ever mean to this query, before the arrangement
        has been taken apart and the question could not be asked again. */

        PluckerPolytope<P> *beam = nullptr;
        if (complex->getPolytopeCount() > 0)
        {
            beam = new PluckerPolytope<P>(*complex->getPolytope(0));
        }

        for (auto occluder = myOccluders.begin(); occluder != myOccluders.end(); occluder++)
        {
            if (VisibilitySolver<P, S>::isTimedOut())
            {
                V_TRACE("FAILURE: time budget over");
                delete beam;
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
                    /* BLOCKED is nought in that enum, which is why it is compared and never just tested. */
                    if (myPolytope != nullptr && isBlockedBySilhouette(myPolytope, myPolyhedron, mySilhouette) == BLOCKED)
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

            /* The faces this query can ever see, and the ones it cannot: proved against the
            beam once their three planes exist, before any cell is cut. An unreachable face
            goes with its three planes out of the arrangement, and out of the ray cast. */
            std::map<const SilhouetteMeshFace *, FaceShadow> shadows;
            collectFaceShadows(myPolyhedron, edges, shadows);
            filterUnreachableFaces(shadows, myPolyhedron, beam, edges);

            std::vector<int> intersectedPolytopes;
            classifyAllPolytopes(edges, complex, myPolyhedron, intersectedPolytopes, mySilhouette);

            treatIntersectedPolytopes(edges, shadows, complex, myPolyhedron, intersectedPolytopes, mySilhouette);

            complex->repack();
        }

        delete beam;
        beam = nullptr;
        /* No conclusion before this: the loops above give up on an over budget query, and the
        work they left undone is not an answer. */
        if (VisibilitySolver<P, S>::isTimedOut())
        {
            V_TRACE("FAILURE: time budget over");
            return FAILURE;
        }

        V_TRACE("End of procedure:" << complex->getPolytopeCount() << " left");
        if (complex->getPolytopeCount() == 0)
        {
             V_TRACE("HIDDEN: done");
       
            return VisibilitySolver<P, S>::conclude(HIDDEN);
        }
         V_TRACE("visible: prepare to extract stabbinglines");
       
        bool validated = false;

        for (int i = 0; i < complex->getPolytopeCount(); i++)
        {
            PluckerPolytope<P> *myPolytope = complex->getPolytope(i);
            V_TRACE(" extract stabbinglines for " << myPolytope << " at " << i);

            /* A cell no silhouette proved blocked is not yet a sight line. The patches only carry the faces the
            silhouette kept, and an occluder that is not closed has faces no silhouette edge bounds: a line can
            cross one of them and still miss every patch. The scene is asked, in full, before it is called visible. */
            if (mOptimisedRepresentativeLine
                && VisibilitySolver<P, S>::mQuery->validateLeafPolytope(myPolytope, myPolyhedron) != 1)
            {
                V_TRACE(myPolytope << " survives the silhouettes but its line is not a sight line");
                continue;
            }

            VisibilitySolver<P, S>::extractStabbingLines(myPolyhedron, myPolytope, mTolerance);

            // the cell survived the whole silhouette set and holds a sight line: trace its discontinuity arcs
            ShadowDiscontinuityTracer<P, S>::traceLeafPolytope(myPolyhedron, myPolytope,
                                                               VisibilitySolver<P, S>::mQuery->getEventRegistry(),
                                                               *VisibilitySolver<P, S>::mQuery->getQueryPolygon(0),
                                                               *VisibilitySolver<P, S>::mQuery->getQueryPolygon(1),
                                                               mTolerance,
                                                               VisibilitySolver<P, S>::mQuery->getShadowDiscontinuityArcs());

            validated = true;
        }

        if (!validated)
        {
            V_TRACE("HIDDEN: no surviving polytope holds a sight line");
            return VisibilitySolver<P, S>::conclude(HIDDEN);
        }

        return VisibilitySolver<P, S>::conclude(VISIBLE);
    }
}
