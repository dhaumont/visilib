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
#include <vector>

#include "geometry_position_type.h"
#include "math_combinatorial.h"
#include "math_plucker_6.h"

namespace visilib
{
    /** @brief Container containing the line coordinates in Plucker space, as well as the facets description of each line. 
 
    The facets description of a line is the list of edges incident to the lines in 3D space.
    In Plucker space, it is list of hyperplanes that meet at the Plucker point of the line.
    For efficiency and precision reasons, the PluckerPolyhedron also store if each line is normalized and its relative position relative to the Plucker quadric

    */

    template<class P>
    class PluckerPolyhedron
    {
    public:
        PluckerPolyhedron();
        ~PluckerPolyhedron();

        /**@brief Return the number of point stored in the polyhedron*/
        size_t getLinesCount() const
        {
            return mLines.size();
        }

        /**@brief Return a point stored in the polyhedron 
 
        @param aNumber: the index of the point
        */
        const P& get(size_t aNumber) const
        {
            return mLines[aNumber];
        }

        /**@brief Return the position of a point with respect to the Plucker Quadric
 
        @param aNumber: the index of the point
        */
        GeometryPositionType getQuadricRelativePosition(size_t aNumber) const
        {
            return mQuadricRelativePositions[aNumber];
        }

        /**@brief Returns if the point is normalized
        @param aNumber: the index of the point
        */
        bool isNormalized(size_t aNumber) const
        {
            return mNormalizations[aNumber];
        }

        /**@brief Add a point to the polyhedron
 
        @param aLine: the point in Plucker space
        @param aPosition: the relative position of the point with respect to the Plucker Quadric
        @param aLine: the normalization status of the point
        @param tolerance: tolerancee used for verification of the data structure

        */
        template<class S> size_t add(const P& aLine, GeometryPositionType aPosition, bool aNormalization, S tolerance);

        /**@brief Return the facet description of a point*/
        std::vector<size_t>& getFacetsDescription(size_t v) 
        { 
            return mFacetsDescription[v];
        }

        /** @brief Add a facet (index of an hyperplane) to the combinatorial description of a vertex
 
        @param aVertice: the index of the vertex
        @param aFace: the index of the facet (index of the hyperplane in mLines)
        */
        void addToFacetsDescription(size_t aVertice, size_t aFace);

        /** @brief Initialize the combinatorial description of a vertex

        @param aVerticeDestination: the index of the vertex
        @param source: the list of facets (index of hyperplane in mLines)
        */
        void initFacetsDescription(size_t aVerticeDestination, const std::vector<size_t> & source);


        /** @brief Check the combinatorial description of a vertex defined at the intersection of an edge and an hyperplane

        The combinatorial description of the edge is computed by performing the intersection of the facet description of its two vertices. 

        @param aVerticeDestination: the index of the new vertex
        @param aVerticeSource1: the first vertex index of the edge
        @param aVerticeSource2: the second vertex index of the edge
        @param aReplacement: the index of the splitting hyperplane
        */
        bool checkFacetsDescription(size_t aVerticeDestination, size_t aVerticeSource1, size_t aVerticeSource2, size_t aReplacement);

        /** @brief Return true if the container contains another line with the same facet description*/
        bool containsOtherLinesWithSameFacetsDescription(size_t duplicatedLine) const;

        void resize(size_t size);

    private:
        /**@brief Check that the vertex has a valid facet description 
 
        The check verify that the facet description is sorted
        */
        bool isValid(size_t aVertice);

        std::vector<P> mLines;                                        /** < @brief The list of Plucker points (Plucker vertices and hyperplanes). */
        std::vector<GeometryPositionType> mQuadricRelativePositions;  /** < @brief The relative position of the Plucker point relative to the Plucker quadric*/
        std::vector<bool> mNormalizations;                            /** < @brief The normalization status of each Plucker point*/
        //std::vector<std::vector<size_t> > mFacetsDescription;          /** < @brief The facet description (list of all the hyperplanes intersecting at that point) of each Plucker point mLines*/
    };

    template<class P>
    template<class S>
    size_t PluckerPolyhedron<P>::add(const P& aLine, GeometryPositionType aPosition, bool aNormalization, S tolerance)
    {
        V_ASSERT(!MathPredicates::isZero(aLine, tolerance));
        V_ASSERT(!aNormalization || MathPredicates::isNormalized(aLine, tolerance));
        mLines.push_back(aLine);
        mQuadricRelativePositions.push_back(aPosition);
        mNormalizations.push_back(aNormalization);
        mFacetsDescription.push_back(std::vector<size_t>());
        return mLines.size() - 1;
    }

    template<class P>
    void PluckerPolyhedron<P>::resize(size_t size)
    {
        mLines.resize(size);
        mQuadricRelativePositions.resize(size);
        mNormalizations.resize(size);
        mFacetsDescription.resize(size);
    }

    template<class P>
    PluckerPolyhedron<P>::PluckerPolyhedron()
    {
    }

    template<class P>
    PluckerPolyhedron<P>::~PluckerPolyhedron()
    {
    }

    template<class P>
    void PluckerPolyhedron<P>::addToFacetsDescription(size_t aVertice, size_t aFace)
    {
        V_ASSERT(isValid(aVertice));

        size_t length = mFacetsDescription[aVertice].size();

        mFacetsDescription[aVertice].resize(length + 1);
        size_t i = length;

        for (; i > 0 && mFacetsDescription[aVertice][i - 1] > aFace; --i)
        {
            mFacetsDescription[aVertice][i] = mFacetsDescription[aVertice][i - 1];
        }
        mFacetsDescription[aVertice][i] = aFace;

        V_ASSERT(isValid(aVertice));
    }

    template<class P>
    void PluckerPolyhedron<P>::initFacetsDescription(size_t aVerticeDestination, const std::vector<size_t> & source)
    {
        //V_ASSERT(mFacetsDescription.size() == aVerticeDestination);

        // a Source is sortded in place, so it can be modified by the function
        std::vector<size_t>& facets = getFacetsDescription(aVerticeDestination);
        for (size_t i = 0; i < source.size(); i++)
        {
            facets.push_back(source[i]);
        }
        std::sort(facets.begin(), facets.end());

        V_ASSERT(isValid(aVerticeDestination));
    }

    template<class P>
    bool PluckerPolyhedron<P>::checkFacetsDescription(size_t aVerticeDestination, size_t aVerticeSource1, size_t aVerticeSource2, size_t aReplacement)
    {
        if (aVerticeDestination >= mFacetsDescription.size())
        {
            return false;
        }

        // sanity checks
        if (!MathCombinatorial::haveAtLeastNCommonFacets(getFacetsDescription(aVerticeSource1), getFacetsDescription(aVerticeSource2)))
            return false;

#ifdef FULL_VALIDITY_CHECK
        if (containsOtherLinesWithSameFacetsDescription(aVerticeDestination))
            return false;
#endif
        if (!isValid(aVerticeDestination))
            return false;
        if (!isValid(aVerticeSource1))
            return false;
        if (!isValid(aVerticeSource2))
            return false;
        if (!MathCombinatorial::haveAtLeastNCommonFacets(getFacetsDescription(aVerticeSource1), getFacetsDescription(aVerticeDestination)))
            return false;
        if (!MathCombinatorial::haveAtLeastNCommonFacets(getFacetsDescription(aVerticeSource2), getFacetsDescription(aVerticeDestination)))
            return false;
        return true;
    }

    template<class P>
    bool PluckerPolyhedron<P>::isValid(size_t aVertice)
    {
        std::vector<size_t>& facets = getFacetsDescription(aVertice);
        return std::is_sorted(facets.begin(), facets.end());
    }
 
    template<class P>
    bool PluckerPolyhedron<P>::containsOtherLinesWithSameFacetsDescription(size_t aVerticeDestination) const
    {
        const std::vector<size_t>& source = mFacetsDescription[aVerticeDestination];
        for (size_t i = 0; i < getLinesCount(); i++)
        {
            if (i != aVerticeDestination)
            {
                const std::vector<size_t>& facets = mFacetsDescription[i];
                if (source.size() == facets.size() && std::equal(source.begin(), source.end(), facets.begin()))
                    return true;
            }
        }
        return false;
    }
}
