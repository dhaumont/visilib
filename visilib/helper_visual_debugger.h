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

#include "math_vector_3.h"
#include "visilib.h"
#include <fstream>

namespace visilib
{
    /** @brief Container for debugging information collected during visibility determination operations

    The container allows to store the stabbing lines, the edges processed and the occluder triangle processed.
    The usage of this container is optional.
    An helperDebugVisualization instance must be attached to the different classes such that they can store debugging information
    */
    class HelperVisualDebugger
    {
    public:
        HelperVisualDebugger();

        ~HelperVisualDebugger();

        /** @brief Clear the container */
        void clear();

        /** @brief Add an occluder edge that has been processed */
        void addRemovedEdge(const MathVector3f& aBegin, const MathVector3f& anEnd);

        /** @brief Add a stabbing line in 3D space */
        void addStabbingLine(const MathVector3f& aBegin, const MathVector3f& anEnd);

        /** @brief Add a sampling line, used to collect occluders, in 3D space */
        void addSamplingLine(const MathVector3f& aBegin, const MathVector3f& anEnd);

        /** @brief Add extreaml stabbing line of a polytope, in 3D space */
        void addExtremalStabbingLine(const MathVector3f& aBegin, const MathVector3f& anEnd);

        /** @brief Add an occluder triangle that has been processed*/
        void addRemovedTriangle(const MathVector3f& a1, const MathVector3f& a2, const MathVector3f& a3);

        /** @brief Return a mesh containing the removed triangles for visualization purposes*/
        const GeometryTriangleMeshDescription& getRemovedTriangles() const;

        /** @brief Return a mesh containing the stabbing lines for visualization purposes*/
        const GeometryTriangleMeshDescription& getStabbingLines() const;

        /** @brief Return a mesh containing the extremal stabbing lines for visualization purposes*/
        const GeometryTriangleMeshDescription& getExtremalStabbingLines() const;

        /** @brief Return a mesh containing the sampling lines for visualization purposes*/
        const GeometryTriangleMeshDescription& getSamplingLines() const;

        /** @brief Return a mesh containing the removed edges for visualization purposes*/
        const GeometryTriangleMeshDescription& getRemovedEdges() const;


#ifdef OUTPUT_DEBUG_FILE
        std::ofstream& getDebugOutput() { return mDebugOutput; };
#endif       
    private:
        std::vector<MathVector3f> mSplittingEdgeArray;   /** <@brief The array of edges that have been processed*/
        std::vector<MathVector3f> mRemovedTriangleArray; /** <@brief The array of triangles that have been processed*/
        std::vector<MathVector3f> mStabbingLineArray;    /** <@brief The array of stabbing lines that have been processed*/
        std::vector<MathVector3f> mSamplingLineArray;   /** <@brief The array of sampling lines that have been processed*/
        std::vector<MathVector3f> mExtremalStabbingLineArray;   /** <@brief The array of extremal stabbbing lines that have been processed*/

        GeometryTriangleMeshDescription mSplittingEdgeInformation;     /** <@brief The mesh containing the removed edges for visualization purposes*/
        GeometryTriangleMeshDescription mRemovedTriangleInformation;   /** <@brief The mesh containing the removed triangles for visualization purposes*/
        GeometryTriangleMeshDescription mStabbingLineInformation;      /** <@brief The mesh containing the removed stabbing lines for visualization purposes*/
        GeometryTriangleMeshDescription mSamplingLineInformation;      /** <@brief The mesh containing the removed sampling lines for visualization purposes*/
        GeometryTriangleMeshDescription mExtremalStabbingLineInformation;   /** <@brief The mesh containing the extremal stabbing lines for visualization purposes*/

#ifdef OUTPUT_DEBUG_FILE
        std::ofstream mDebugOutput;
#endif
    };

    inline HelperVisualDebugger::HelperVisualDebugger()
    {

#ifdef OUTPUT_DEBUG_FILE
        mDebugOutput.open("debug.txt");
#endif
        clear();
    }

    inline HelperVisualDebugger::~HelperVisualDebugger()
    {
#ifdef OUTPUT_DEBUG_FILE
        mDebugOutput.close();
#endif
    }

    inline void HelperVisualDebugger::addRemovedEdge(const MathVector3f& aBegin, const MathVector3f& anEnd)
    {
        mSplittingEdgeArray.push_back(aBegin);
        mSplittingEdgeArray.push_back(anEnd);

        mSplittingEdgeInformation.vertexCount = mSplittingEdgeArray.size();
        mSplittingEdgeInformation.vertexArray = &mSplittingEdgeArray[0];
    }

    inline void HelperVisualDebugger::addStabbingLine(const MathVector3f& aBegin, const MathVector3f& anEnd)
    {
        mStabbingLineArray.push_back(aBegin);
        mStabbingLineArray.push_back(anEnd);

        mStabbingLineInformation.vertexCount = mStabbingLineArray.size();
        mStabbingLineInformation.vertexArray = &mStabbingLineArray[0];
    }

    inline void HelperVisualDebugger::addSamplingLine(const MathVector3f& aBegin, const MathVector3f& anEnd)
    {
        mSamplingLineArray.push_back(aBegin);
        mSamplingLineArray.push_back(anEnd);

        mSamplingLineInformation.vertexCount = mSamplingLineArray.size();
        mSamplingLineInformation.vertexArray = &mSamplingLineArray[0];
    }

    inline void HelperVisualDebugger::addExtremalStabbingLine(const MathVector3f& aBegin, const MathVector3f& anEnd)
    {
        mExtremalStabbingLineArray.push_back(aBegin);
        mExtremalStabbingLineArray.push_back(anEnd);

        mExtremalStabbingLineInformation.vertexCount = mExtremalStabbingLineArray.size();
        mExtremalStabbingLineInformation.vertexArray = &mExtremalStabbingLineArray[0];
    }
    inline void  HelperVisualDebugger::addRemovedTriangle(const MathVector3f& a1, const MathVector3f& a2, const MathVector3f& a3)
    {
        mRemovedTriangleArray.push_back(a1);
        mRemovedTriangleArray.push_back(a2);
        mRemovedTriangleArray.push_back(a3);

        mRemovedTriangleInformation.vertexCount = mRemovedTriangleArray.size();
        mRemovedTriangleInformation.vertexArray = &mRemovedTriangleArray[0];
    }

    inline void HelperVisualDebugger::clear()
    {

#ifdef OUTPUT_DEBUG_FILE
        mDebugOutput.close();
        mDebugOutput.open("debug.txt");
        mDebugOutput.clear();
#endif
        mSplittingEdgeArray.clear();
        mRemovedTriangleArray.clear();
        mStabbingLineArray.clear();
        mSamplingLineArray.clear();
        mExtremalStabbingLineArray.clear();

        mSplittingEdgeInformation.vertexCount = 0;
        mSplittingEdgeInformation.vertexArray = nullptr;

        mSamplingLineInformation.vertexCount = 0;
        mSamplingLineInformation.vertexArray = nullptr;

        mStabbingLineInformation.vertexCount = 0;
        mStabbingLineInformation.vertexArray = nullptr;

        mRemovedTriangleInformation.vertexCount = 0;
        mRemovedTriangleInformation.vertexArray = nullptr;

        mExtremalStabbingLineInformation.vertexCount = 0;
        mExtremalStabbingLineInformation.vertexArray = nullptr;

    }

    inline const GeometryTriangleMeshDescription& HelperVisualDebugger::getRemovedTriangles() const
    {
        return mRemovedTriangleInformation;
    }

    inline const GeometryTriangleMeshDescription& HelperVisualDebugger::getStabbingLines() const
    {
        return mStabbingLineInformation;
    }

    inline const GeometryTriangleMeshDescription& HelperVisualDebugger::getExtremalStabbingLines() const
    {
        return mExtremalStabbingLineInformation;
    }

    inline const GeometryTriangleMeshDescription& HelperVisualDebugger::getSamplingLines() const
    {
        return mSamplingLineInformation;
    }

    inline const GeometryTriangleMeshDescription& HelperVisualDebugger::getRemovedEdges() const
    {
        return mSplittingEdgeInformation;
    }
}