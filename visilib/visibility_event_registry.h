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

#include <array>
#include <map>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "visibility_event.h"

namespace visilib
{
    /** @brief The provenance of the hyperplanes of the polyhedron: which 3D edge each of them carries.

    The polyhedron knows its hyperplanes only as Plucker points; the events need the edges behind them, and need to know
    when two hyperplanes carry the same physical edge (a mesh edge shared by two silhouettes, say). The registry gives
    each distinct edge a feature id, the same for both orientations of the edge, and maps the hyperplane indices to it.
    A hyperplane the registry does not know (a vertex of a cell, or an edge registered by no one) has no feature.

    Distinct edges can lie on the same line: the sides of a subdivided polygon, cut in collinear pieces. Their hyperplanes
    are the same plane of the Plucker space, and a cell bounded by one of them is bounded by all: the critical lines that
    plane carries may lean on any of the pieces, not only on the one whose hyperplane the cell happens to record. The
    registry gathers the occluder edges of one line in a group, so that an event can be read on every piece of it.*/
    class VisibilityEventRegistry
    {
    public:
        static constexpr size_t cNoFeature = static_cast<size_t>(-1);

        void clear()
        {
            mFeatures.clear();
            mHyperplanes.clear();
            mFeatureOfHyperplane.clear();
            mFeatureOfEdge.clear();
            mGroupOfFeature.clear();
            mGroups.clear();
        }

        /** @brief Record that the hyperplane anIndex carries aFeature. Return the feature id.*/
        size_t registerHyperplane(size_t anIndex, const VisibilityFeature& aFeature)
        {
            const size_t id = getOrCreateFeature(aFeature);
            if (mFeatureOfHyperplane.insert(std::make_pair(anIndex, id)).second)
            {
                mHyperplanes.push_back(anIndex);
            }
            return id;
        }

        /** @brief The feature id of the hyperplane anIndex, cNoFeature if it is not registered.*/
        size_t getFeatureId(size_t anIndex) const
        {
            auto iter = mFeatureOfHyperplane.find(anIndex);
            return iter == mFeatureOfHyperplane.end() ? cNoFeature : iter->second;
        }

        bool hasHyperplane(size_t anIndex) const
        {
            return mFeatureOfHyperplane.find(anIndex) != mFeatureOfHyperplane.end();
        }

        const VisibilityFeature& getFeature(size_t aFeatureId) const
        {
            return mFeatures[aFeatureId];
        }

        /** @brief The collinear group of a feature: every registered edge on its line (itself included), by feature id.*/
        const std::vector<size_t>& getCollinearFeatures(size_t aFeatureId) const
        {
            return mGroups[mGroupOfFeature[aFeatureId]];
        }

        /** @brief The id of the collinear group of a feature: two features with the same group lie on the same line.*/
        size_t getGroupId(size_t aFeatureId) const
        {
            return mGroupOfFeature[aFeatureId];
        }

        /** @brief The indices of the registered hyperplanes, in the order of their registration.*/
        const std::vector<size_t>& getHyperplanes() const
        {
            return mHyperplanes;
        }

        size_t getFeatureCount() const
        {
            return mFeatures.size();
        }

    private:
        typedef std::array<double, 7> EdgeKey; // kind, then the two endpoints in lexicographic order

        static EdgeKey makeKey(const VisibilityFeature& aFeature)
        {
            const MathVector3d& a = aFeature.mEndpoints[0];
            const MathVector3d& b = aFeature.mEndpoints[1];
            const bool swap = std::make_tuple(b.x, b.y, b.z) < std::make_tuple(a.x, a.y, a.z);
            const MathVector3d& first = swap ? b : a;
            const MathVector3d& second = swap ? a : b;
            return EdgeKey{ double(aFeature.mKind), first.x, first.y, first.z, second.x, second.y, second.z };
        }

        size_t getOrCreateFeature(const VisibilityFeature& aFeature)
        {
            const EdgeKey key = makeKey(aFeature);
            auto iter = mFeatureOfEdge.find(key);
            if (iter != mFeatureOfEdge.end())
            {
                return iter->second;
            }
            const size_t id = mFeatures.size();
            mFeatures.push_back(aFeature);
            mFeatureOfEdge[key] = id;

            // the group of the edges on its line: only the occluder edges are gathered, those of the sources stay alone
            size_t group = mGroups.size();
            if (aFeature.mKind == VisibilityFeature::OCCLUDER_EDGE)
            {
                for (size_t other = 0; other < id; other++)
                {
                    if (mFeatures[other].mKind == VisibilityFeature::OCCLUDER_EDGE && isCollinear(mFeatures[other], aFeature))
                    {
                        group = mGroupOfFeature[other];
                        break;
                    }
                }
            }
            if (group == mGroups.size())
            {
                mGroups.push_back(std::vector<size_t>());
            }
            mGroups[group].push_back(id);
            mGroupOfFeature.push_back(group);
            return id;
        }

        /** @brief Whether both endpoints of b lie on the line of a, within a relative tolerance of the length of a.*/
        static bool isCollinear(const VisibilityFeature& a, const VisibilityFeature& b)
        {
            const MathVector3d d = a.mEndpoints[1] - a.mEndpoints[0];
            const double length = d.getNorm();
            if (!(length > 0.0))
            {
                return false;
            }
            for (int i = 0; i < 2; i++)
            {
                const double distance = MathVector3d::cross(d, b.mEndpoints[i] - a.mEndpoints[0]).getNorm() / length;
                if (distance > cCollinearTolerance * length)
                {
                    return false;
                }
            }
            return true;
        }

        static constexpr double cCollinearTolerance = 1e-6; /**< @brief The distance, relative to the length, under which a point is on a line*/

        std::vector<VisibilityFeature> mFeatures;                  /**< @brief The distinct edges, by feature id*/
        std::vector<size_t> mHyperplanes;                          /**< @brief The registered hyperplane indices*/
        std::unordered_map<size_t, size_t> mFeatureOfHyperplane;   /**< @brief The feature id of each registered hyperplane*/
        std::map<EdgeKey, size_t> mFeatureOfEdge;                  /**< @brief The feature id of each distinct edge*/
        std::vector<size_t> mGroupOfFeature;                       /**< @brief The collinear group of each feature*/
        std::vector<std::vector<size_t>> mGroups;                  /**< @brief The features of each collinear group*/
    };
}
