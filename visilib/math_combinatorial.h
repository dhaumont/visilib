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

#include <vector>
#include <algorithm>
#include "visilib_core.h"

namespace visilib
{
    /** @brief Provides functions operating on the facets description of Plucker points.*/

    class MathCombinatorial
    {
    public:
        static bool haveAtLeastNCommonFacets(const std::vector<size_t>& aFacetsDescription1, const std::vector<size_t>& aFacetsDescription2, size_t n = 3);
        static void initFacets(const std::vector<size_t>& aFacetsDescription1, const std::vector<size_t>& aFacetsDescription2, size_t anHyperplane, std::vector<size_t>& result);
        static void initFacets(const std::vector<size_t>& aFacetsDescription1, const std::vector<size_t>& aFacetsDescription2, std::vector<size_t>& result);
        static bool hasFacet(const std::vector<size_t>& facets, size_t aFace);
    };

    /** @brief Determine if the intersection of the facets lists of the two Plucker Points have at least n common elements.

    This function is used to determine if two vertices of a Polytope have to be linked by an edge during polytope manipulation routines.
    n is 3 by default for Plucker line in 3D space.*/

    inline bool MathCombinatorial::haveAtLeastNCommonFacets(const std::vector<size_t>& aFacetsDescription1, const std::vector<size_t>& aFacetsDescription2, size_t n)
    {
        V_ASSERT(std::is_sorted(aFacetsDescription1.begin(), aFacetsDescription1.end()));
        V_ASSERT(std::is_sorted(aFacetsDescription2.begin(), aFacetsDescription2.end()));

        std::vector<size_t> result(aFacetsDescription1.size());
        auto iter = std::set_intersection(aFacetsDescription1.begin(), aFacetsDescription1.end(), aFacetsDescription2.begin(), aFacetsDescription2.end(), result.begin());
        size_t size = iter - result.begin();

        return size >= n;
    }

    /** @brief Initialize a new facets description result from the two input facets description aFacetsDescription1 an aFacetsDescription1 of two given vertices and an additional
    facet anHyperplane.
    The resulting facets description is a sorted list of facets, containing the intersection set of the two input facets descriptions list and the additional facet anHyperplane.
    Remark: the precondition is that the input facets list are sorted, that they have at least n common facets and that anHyperplane is greater than the facets of the lists.*/

    inline void MathCombinatorial::initFacets(const std::vector<size_t>& aFacetsDescription1, const std::vector<size_t>& aFacetsDescription2, size_t anHyperplane, std::vector<size_t>& result)
    {
        V_ASSERT(MathCombinatorial::haveAtLeastNCommonFacets(aFacetsDescription1, aFacetsDescription2));
        bool requireSorting = anHyperplane <= aFacetsDescription1[aFacetsDescription1.size() - 1]
         || anHyperplane <= aFacetsDescription2[aFacetsDescription2.size() - 1];

        result.resize(aFacetsDescription1.size());
        auto iter = std::set_intersection(aFacetsDescription1.begin(), aFacetsDescription1.end(), aFacetsDescription2.begin(), aFacetsDescription2.end(), result.begin());
        size_t size = iter - result.begin();
        result.resize(size + 1);
        result[size] = anHyperplane;
        if (requireSorting)
        {
            std::sort(result.begin(), result.end());
        }
        V_ASSERT(std::is_sorted(result.begin(), result.end()));
    }

    /** @brief Initialize a new facets description result from the two input facets description aFacetsDescription1 an aFacetsDescription1 of two given vertices
    The resulting facets description is a sorted list of facets, containing the intersection set of the two input facets descriptions list and the additional facet anHyperplane.
    Remark: the precondition is that the input facets list are sorted .*/

    inline void MathCombinatorial::initFacets(const std::vector<size_t>& aFacetsDescription1, const std::vector<size_t>& aFacetsDescription2, std::vector<size_t>& result)
    {
        int i = 0, j = 0;
        int n1 = aFacetsDescription1.size();
        int n2 = aFacetsDescription2.size();

        while (i < n1 && j < n2)
        {
            if (aFacetsDescription1[i] < aFacetsDescription2[j])
            {
                aFacetsDescription1[i++];
            }
            else if (aFacetsDescription1[i] > aFacetsDescription2[j])
            {
                result.push_back(aFacetsDescription2[j++]);
            }
            else
            {
                i++;j++;
            }
        }

        while (i < n1)
        {
            result.push_back(aFacetsDescription1[i++]);
        }
        while (j < n2)
        {
            result.push_back(aFacetsDescription2[j++]);
        }
    }

    /** @brief Determine if the sorted facets description contains aFace. Behaviour is only defined if aFacetsDescription is sorted.*/

    inline bool MathCombinatorial::hasFacet(const std::vector<size_t> & aFacetsDescription, size_t aFace)
    {
        V_ASSERT(std::is_sorted(aFacetsDescription.begin(), aFacetsDescription.end()));

        return std::binary_search(aFacetsDescription.begin(), aFacetsDescription.end(), aFace);;
    }
}
