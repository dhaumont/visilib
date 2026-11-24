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

#include <map>
#include <set>
#include <utility>
#include <vector>

#include "geometry_convex_polygon.h"
#include "math_combinatorial.h"
#include "math_conic_2.h"
#include "math_frame_3.h"
#include "math_geometry.h"
#include "math_plucker_6.h"
#include "math_plucker_traits.h"
#include "plucker_polyhedron.h"
#include "plucker_polytope.h"
#include "shadow_event_family.h"
#include "visibility_event.h"
#include "visibility_event_registry.h"

namespace visilib
{
    /** @brief Collects the shadow discontinuities of the visible cells of a query, and traces them on the receiver plane.

    Query polygon 0 is the light, query polygon 1 the receiver. Each visible cell of the occlusion tree is a set of free
    lines from the light to the receiver, bounded by the hyperplanes of edges. A visibility event of the cell is a set of
    its hyperplanes whose common lines form a one parameter family: three of them in the general space, the lines
    meeting three edges (cCommonFacetsForEdge), one in the directional family, the rays meeting one edge. The free part
    of the event is the range of the family the cell holds: inside every other facet, and meeting the real segments of
    the edges, not only their lines.

    The candidates are the subsets of cCommonFacetsForEdge facets of each vertex of the cell: a face of the cell carrying
    an event has vertices, and they lie on its facets. Each candidate is read through the registry of the edges:
    - an edge of the receiver: the lines run along the boundary of the receiver, not a shadow. Left out.
    - two edges sharing a real vertex v: the lines through v meeting the remaining edge (a VE event). The other lines
      meeting both edges lie in their plane, the plane of the light or of a face of an occluder, and are left out.
    - three edges pairwise skew: their transversals (an EEE event).
    - two coplanar edges with no common vertex: the lines meeting both lie in their plane or go through the point where
      they cross, outside one of them. Left out.

    The same event is found by every cell along it, possibly through different hyperplanes of the same edge: the free
    ranges are gathered by the edges of the event (the feature ids of the registry) and joined, so that each maximal
    range gives one arc.*/
    template<class P, class S>
    class ShadowDiscontinuityExtractor
    {
    public:
        /** @brief Prepare the extraction for aLight and aReceiver, the hyperplanes being described by aRegistry.*/
        ShadowDiscontinuityExtractor(const VisibilityEventRegistry* aRegistry, const GeometryConvexPolygon& aLight, const GeometryConvexPolygon& aReceiver)
            : mRegistry(aRegistry),
            mVisibleCellCount(0),
            mCandidateCount(0),
            mDegenerateCount(0)
        {
            const MathVector3d lightCenter = MathGeometry::getGravityCenter(aLight);
            const MathVector3d receiverCenter = MathGeometry::getGravityCenter(aReceiver);
            mOrientation = receiverCenter - lightCenter;

            MathVector3d axis = MathVector3d::X();
            if (aReceiver.getVertexCount() >= 2)
            {
                axis = aReceiver.getVertex(1) - aReceiver.getVertex(0);
            }
            MathVector3d normal = aReceiver.getPlane().mNormal;
            if (MathVector3d::cross(normal, axis).getSquaredNorm() <= 1e-24 * normal.getSquaredNorm() * axis.getSquaredNorm())
            {
                MathVector3d v;
                MathGeometry::getTangentBasis(normal.getNormalized(), axis, v);
            }
            mFrame = MathFrame3d(receiverCenter, normal, axis);
        }

        /** @brief Gather the free events of a visible cell. The cell must be alive: its vertices and facets are read in place.*/
        void addVisibleCell(PluckerPolyhedron<P>* aPolyhedron, PluckerPolytope<P>* aPolytope, const S& aTolerance)
        {
            mVisibleCellCount++;

            /* The cell is bounded by every edge it lies behind, not only by the hyperplanes tight on its vertices: those
            are a part of the facets of the hull of the vertices, and lines that satisfy them alone can leave the cell
            (across an edge of the receiver, say). Every registered edge the cell lies behind is a bound, among them all
            the splits from the root to the cell. */
            std::vector<typename PluckerPolytope<P>::CellFacet> cellFacets;
            bool flat = false;
            aPolytope->getCellFacets(aPolyhedron, aTolerance, mRegistry->getHyperplanes(), cellFacets, flat);

            std::vector<ShadowCellFacet> facets;
            facets.reserve(cellFacets.size());
            for (const auto& f : cellFacets)
            {
                const size_t id = mRegistry->getFeatureId(f.index);
                facets.push_back(ShadowCellFacet{ f.direction, f.location, f.sign, id == VisibilityEventRegistry::cNoFeature ? VisibilityEventRegistry::cNoFeature : mRegistry->getGroupId(id) });
            }

            // the directional family: its rays all have the direction of any vertex of the cell
            MathVector3d rayDirection = MathVector3d::Zero();
            if constexpr (MathPluckerCellTraits<P>::cAllLinesReal)
            {
                if (!aPolytope->getVertices().empty())
                {
                    rayDirection = convert<MathVector3d>(aPolyhedron->get(*aPolytope->getVertices().begin()).getDirection());
                    if (rayDirection.dot(mOrientation) < 0.0)
                    {
                        rayDirection *= -1.0;
                    }
                }
            }

            std::set<std::vector<size_t>> candidates;
            for (size_t vertex : aPolytope->getVertices())
            {
                MathCombinatorial::forEachSubset(aPolyhedron->getFacetsDescription(vertex), MathPluckerCellTraits<P>::cCommonFacetsForEdge,
                    [&](const std::vector<size_t>& subset) { candidates.insert(subset); });
            }

            for (const std::vector<size_t>& candidate : candidates)
            {
                mCandidateCount++;
                std::vector<size_t> features;
                bool known = true;
                for (size_t hyperplane : candidate)
                {
                    const size_t id = mRegistry->getFeatureId(hyperplane);
                    if (id == VisibilityEventRegistry::cNoFeature)
                    {
                        known = false;
                        break;
                    }
                    features.push_back(id);
                }
                if (!known)
                {
                    continue;
                }

                if constexpr (MathPluckerCellTraits<P>::cAllLinesReal)
                {
                    const size_t group = mRegistry->getGroupId(features[0]);
                    std::vector<ShadowCellFacet> others;
                    for (const ShadowCellFacet& facet : facets)
                    {
                        if (facet.mGroup != group)
                        {
                            others.push_back(facet);
                        }
                    }
                    for (size_t piece : mRegistry->getCollinearFeatures(features[0]))
                    {
                        addDirectionalEvent(piece, rayDirection, others);
                    }
                }
                else
                {
                    addEvent(features, facets);
                }
            }
        }

        /** @brief Trace the gathered events on the receiver plane: append one arc per maximal free range.*/
        void getArcs(std::vector<ShadowDiscontinuityArc>& arcs) const
        {
            for (const auto& entry : mEvents)
            {
                const EventRecord& record = entry.second;
                std::vector<std::pair<double, double>> ranges = record.mRanges;
                std::sort(ranges.begin(), ranges.end());

                std::vector<std::pair<double, double>> joined;
                for (const auto& range : ranges)
                {
                    if (!joined.empty() && range.first <= joined.back().second + 1e-9)
                    {
                        joined.back().second = std::max(joined.back().second, range.second);
                    }
                    else
                    {
                        joined.push_back(range);
                    }
                }

                for (const auto& range : joined)
                {
                    ShadowDiscontinuityArc arc;
                    if (traceArc(record, range.first, range.second, arc))
                    {
                        arcs.push_back(arc);
                    }
                }
            }
        }

        const MathFrame3d& getFrame() const { return mFrame; }
        size_t getVisibleCellCount() const { return mVisibleCellCount; }
        size_t getCandidateCount() const { return mCandidateCount; }
        size_t getDegenerateCount() const { return mDegenerateCount; }
        size_t getEventCount() const { return mEvents.size(); }

    private:
        /** @brief An event and the free ranges of its family gathered over the cells.*/
        struct EventRecord
        {
            ShadowDiscontinuityArc::Event mEvent;
            ShadowEventFamily mFamily;
            std::vector<MathVector3d> mEdges;   // the three lines (two points each) the conic of an EEE event is made of, or the edges of a VE event
            std::vector<std::pair<double, double>> mRanges;
        };

        /** @brief The key of an event: its kind, then the ids and coordinates that make its family.*/
        typedef std::vector<double> EventKey;

        /** @brief Read the event of three facets on every piece of their lines.

        A facet carries a line, and the edges of that line (its collinear group) are all pieces of the boundary the
        facet stands for: the event is read for each choice of one piece per facet, the clipping keeping the lines that
        lean on the piece itself. Two facets on the same line are one constraint, not two: their lines meet along all
        of it, and the triple is not an event.*/
        void addEvent(const std::vector<size_t>& someFeatures, const std::vector<ShadowCellFacet>& someFacets)
        {
            size_t group[3];
            for (int i = 0; i < 3; i++)
            {
                if (mRegistry->getFeature(someFeatures[i]).mFlat)
                {
                    return; // an edge inside a flat region of an occluder: the lines leaning on it cross the occluder
                }
                group[i] = mRegistry->getGroupId(someFeatures[i]);
            }
            if (group[0] == group[1] || group[0] == group[2] || group[1] == group[2])
            {
                return;
            }

            /* The lines of the event lie on the hyperplanes of its own lines: on the three facets, and on the hyperplanes
            of the other pieces of the same lines, which only the rounding of the coordinates keeps from being the same
            planes. Clipping the family by them would read that rounding as a side: they are left out. */
            std::vector<ShadowCellFacet> facets;
            facets.reserve(someFacets.size());
            for (const ShadowCellFacet& facet : someFacets)
            {
                if (facet.mGroup != group[0] && facet.mGroup != group[1] && facet.mGroup != group[2])
                {
                    facets.push_back(facet);
                }
            }

            for (size_t a : mRegistry->getCollinearFeatures(someFeatures[0]))
            {
                for (size_t b : mRegistry->getCollinearFeatures(someFeatures[1]))
                {
                    for (size_t c : mRegistry->getCollinearFeatures(someFeatures[2]))
                    {
                        // the event of the pieces, identified by their ids in a canonical order
                        std::vector<size_t> pieces = { a, b, c };
                        std::sort(pieces.begin(), pieces.end());
                        addSegmentEvent(pieces, facets);
                    }
                }
            }
        }

        void addSegmentEvent(const std::vector<size_t>& someFeatures, const std::vector<ShadowCellFacet>& someFacets)
        {
            const VisibilityFeature* f[3];
            for (int i = 0; i < 3; i++)
            {
                f[i] = &mRegistry->getFeature(someFeatures[i]);
                if (f[i]->mKind == VisibilityFeature::RECEIVER_EDGE)
                {
                    return;
                }
            }

            // the vertices shared by two of the edges: each one gives the pencil over the third edge
            bool hasVertex = false;
            for (int i = 0; i < 3; i++)
            {
                for (int j = i + 1; j < 3; j++)
                {
                    const int shared = f[i]->getSharedRealEndpoint(*f[j]);
                    if (shared < 0)
                    {
                        continue;
                    }
                    hasVertex = true;
                    const int k = 3 - i - j;
                    const MathVector3d& vertex = f[i]->mEndpoints[shared];
                    const VisibilityFeature& edge = *f[k];
                    if (edge.mSynthetic || edge.getSharedRealEndpoint(*f[i]) >= 0 || edge.getSharedRealEndpoint(*f[j]) >= 0)
                    {
                        // the third edge is padding, or touches the vertex too: the lines lie in a plane of the light or of an occluder
                        continue;
                    }
                    const bool lightVertex = f[i]->mKind == VisibilityFeature::LIGHT_EDGE && f[j]->mKind == VisibilityFeature::LIGHT_EDGE;
                    if (lightVertex && edge.mKind == VisibilityFeature::LIGHT_EDGE)
                    {
                        continue;
                    }
                    if (!lightVertex && (f[i]->mKind == VisibilityFeature::LIGHT_EDGE || f[j]->mKind == VisibilityFeature::LIGHT_EDGE))
                    {
                        // a light edge and an occluder edge sharing a point: the light touches the occluder
                        mDegenerateCount++;
                        continue;
                    }

                    EventKey key = { double(lightVertex ? ShadowDiscontinuityArc::VE_LIGHT_VERTEX : ShadowDiscontinuityArc::VE_OCCLUDER_VERTEX),
                                     vertex.x, vertex.y, vertex.z, double(someFeatures[k]) };
                    auto iter = mEvents.find(key);
                    if (iter == mEvents.end())
                    {
                        EventRecord record;
                        record.mEvent = lightVertex ? ShadowDiscontinuityArc::VE_LIGHT_VERTEX : ShadowDiscontinuityArc::VE_OCCLUDER_VERTEX;
                        if (!record.mFamily.buildVertexPencil(vertex, edge.mEndpoints[0], edge.mEndpoints[1], mOrientation))
                        {
                            mDegenerateCount++;
                            continue;
                        }
                        record.mEdges = { f[i]->mEndpoints[0], f[i]->mEndpoints[1], f[j]->mEndpoints[0], f[j]->mEndpoints[1], edge.mEndpoints[0], edge.mEndpoints[1] };
                        iter = mEvents.insert(std::make_pair(key, record)).first;
                    }
                    iter->second.mFamily.getFreeIntervals(someFacets, cTolerance, iter->second.mRanges);
                }
            }
            if (hasVertex)
            {
                return;
            }

            // no shared vertex: three edges, at most one of them on the light
            size_t lightEdges = 0;
            for (int i = 0; i < 3; i++)
            {
                if (f[i]->mSynthetic)
                {
                    return;
                }
                lightEdges += f[i]->mKind == VisibilityFeature::LIGHT_EDGE ? 1 : 0;
            }
            if (lightEdges > 1)
            {
                return; // two edges of the light: the lines meeting both lie in its plane
            }

            MathPlucker6d lines[3];
            for (int i = 0; i < 3; i++)
            {
                lines[i] = MathPlucker6d(f[i]->mEndpoints[0], f[i]->mEndpoints[1]);
            }
            for (int i = 0; i < 3; i++)
            {
                for (int j = i + 1; j < 3; j++)
                {
                    if (isCoplanar(lines[i], lines[j]))
                    {
                        mDegenerateCount++;
                        return;
                    }
                }
            }

            EventKey key = { double(ShadowDiscontinuityArc::EEE), double(someFeatures[0]), double(someFeatures[1]), double(someFeatures[2]) };
            auto iter = mEvents.find(key);
            if (iter == mEvents.end())
            {
                EventRecord record;
                record.mEvent = ShadowDiscontinuityArc::EEE;
                const MathVector3d b[3] = { f[0]->mEndpoints[0], f[1]->mEndpoints[0], f[2]->mEndpoints[0] };
                const MathVector3d e[3] = { f[0]->mEndpoints[1], f[1]->mEndpoints[1], f[2]->mEndpoints[1] };
                if (!record.mFamily.buildTransversals(b, e, mOrientation))
                {
                    mDegenerateCount++;
                    return;
                }
                record.mEdges = { b[0], e[0], b[1], e[1], b[2], e[2] };
                iter = mEvents.insert(std::make_pair(key, record)).first;
            }
            iter->second.mFamily.getFreeIntervals(someFacets, cTolerance, iter->second.mRanges);
        }

        void addDirectionalEvent(size_t aFeature, const MathVector3d& aRayDirection, const std::vector<ShadowCellFacet>& someFacets)
        {
            const VisibilityFeature& feature = mRegistry->getFeature(aFeature);
            if (feature.mKind != VisibilityFeature::OCCLUDER_EDGE || feature.mFlat)
            {
                return; // the rays along the outline of the light or of the receiver
            }
            EventKey key = { double(ShadowDiscontinuityArc::DIRECTIONAL), double(aFeature) };
            auto iter = mEvents.find(key);
            if (iter == mEvents.end())
            {
                EventRecord record;
                record.mEvent = ShadowDiscontinuityArc::DIRECTIONAL;
                if (!record.mFamily.buildDirectional(aRayDirection, feature.mEndpoints[0], feature.mEndpoints[1]))
                {
                    mDegenerateCount++;
                    return;
                }
                record.mEdges = { feature.mEndpoints[0], feature.mEndpoints[1] };
                iter = mEvents.insert(std::make_pair(key, record)).first;
            }
            iter->second.mFamily.getFreeIntervals(someFacets, cTolerance, iter->second.mRanges);
        }

        static bool isCoplanar(const MathPlucker6d& a, const MathPlucker6d& b)
        {
            // the Klein form vanishes for lines meeting or parallel, relative to the length of the directions and the size of the moments
            const double scale = a.getDirection().getNorm() * b.getLocation().getNorm() + a.getLocation().getNorm() * b.getDirection().getNorm()
                + a.getDirection().getNorm() * b.getDirection().getNorm();
            return std::fabs(a.dot(b)) <= 1e-12 * scale;
        }

        bool traceArc(const EventRecord& aRecord, double t0, double t1, ShadowDiscontinuityArc& arc) const
        {
            arc.mEvent = aRecord.mEvent;
            arc.mType = aRecord.mEvent == ShadowDiscontinuityArc::EEE ? ShadowDiscontinuityArc::CONIC : ShadowDiscontinuityArc::SEGMENT;
            arc.mFrame = mFrame;
            arc.mT0 = t0;
            arc.mT1 = t1;
            arc.mDirection = aRecord.mFamily.getDirectionCoefficients();
            arc.mLocation = aRecord.mFamily.getLocationCoefficients();
            arc.mEventEdges = aRecord.mEdges;

            if (!arc.getPoint(t0, arc.mEndpoints[0]) || !arc.getPoint(t1, arc.mEndpoints[1]))
            {
                return false;
            }
            if (t1 - t0 < cMinimumRange || (arc.mEndpoints[0] - arc.mEndpoints[1]).getNorm() <= 1e-12 * (1.0 + arc.mEndpoints[0].getNorm()))
            {
                return false; // the family shrinks to a single line on the receiver: a point, not a curve
            }

            arc.mLineCoefficients = MathVector2d::getHomogeneousLine(arc.mEndpoints[0], arc.mEndpoints[1]);
            if (arc.mType == ShadowDiscontinuityArc::CONIC)
            {
                MathPlucker6d lines[3];
                for (int i = 0; i < 3; i++)
                {
                    lines[i] = MathPlucker6d(aRecord.mEdges[2 * i], aRecord.mEdges[2 * i + 1]);
                }
                arc.mConic = MathConic2d::fromTransversalsOfThreeLines(lines, mFrame);
            }
            return true;
        }

        static constexpr double cMinimumRange = 1e-6; /**< @brief The range of the parameter, along an edge of unit length, under which an arc is a point*/
        static constexpr double cTolerance = 1e-9; /**< @brief How far outside a facet, as a distance, a line is still taken as in*/

        const VisibilityEventRegistry* mRegistry;
        MathFrame3d mFrame;
        MathVector3d mOrientation;
        std::map<EventKey, EventRecord> mEvents;
        size_t mVisibleCellCount;
        size_t mCandidateCount;
        size_t mDegenerateCount;
    };
}
