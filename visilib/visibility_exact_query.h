/*
Visilib, an open source library for exact visibility computation.
Copyright(C) 2020 by Denis Haumont

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

#include <unordered_map>
#include <unordered_set>

#include "geometry_convex_polygon.h"
#include "silhouette_mesh_face.h"
#include "geometry_occluder_set.h"
#include "helper_debug_visualisation.h"
#include "helper_statistic_collector.h"
#include "math_geometry.h"
#include "math_predicates.h"
#include "math_plucker_2.h"
#include "math_vector_3.h"
#include "plucker_polytope.h"
#include "plucker_polytope_builder.h"
#include "plucker_polytope_complex.h"
#include "visibility_aperture_finder.h"
#include "silhouette_container.h"
#include "silhouette_processor.h"
#include "visilib.h"
#include "visilib_core.h"
#include "silhouette_container_embree.h"

namespace visilib
{
    class HelperDebugVisualisation;
    class HelperStatisticCollector;

    class GeometryConvexPolygon;
    class SilhouetteMeshFace;
    class GeometryOccluderSet;
    struct Ray;
    class VisibilitySilhouette;
    class SilhouetteProcessor;

    template <class S>
    class PluckerPolytope;

    /**@brief Exact visibility query interface class

    This interface is the protoype definition of an exact visibility query between two source polygons
    */
    class IVisibilityExactQuery
    {
    public:
        virtual ~IVisibilityExactQuery() {};

        virtual void attachVisualisationDebugger(VisualDebugger* aDebugger) = 0;
        virtual VisibilityResult arePolygonsVisible(const float* vertices0, size_t numVertices0, const float* vertices1, size_t numVertices1) = 0;
        virtual void displayStatistic() = 0;
    };


    /**@brief Exact Visibility solver implementation class

    This class manages the computation of exact visibility query between two polygons through a set of occluders contained in a scene.
    */
    template<class P, class S>
    class VisibilityExactQuery : public IVisibilityExactQuery
    {
    public:

        VisibilityExactQuery(GeometryOccluderSet* aScene, const QueryConfiguration& aConfiguration, S tolerance);

        virtual ~VisibilityExactQuery();

        /**@brief Return the source polygons at the given index */
        GeometryConvexPolygon* getQueryPolygon(size_t i)
        {
            return mQueryPolygon[i];
        }

        /**@brief Return the polytope complex encoding the occlusion tree */
        PluckerPolytopeComplex<P>* getComplex()
        {
            return mComplex;
        }

        /**@brief Performs the intersection between a segment and the geometry of the scene.

        When an intersection is found, the silhouette are extracted from the intersected face and the silhouette is linked to the polytope
        The links is used for occluder selection.
        */
        bool findSceneIntersection(const MathVector3d& aBegin, const MathVector3d& anEnd, std::set<SilhouetteMeshFace*>& intersectedFaces, PluckerPolytope<P>* aPolytope = nullptr);

        void extractAllSilhouettes();

        /**@brief Given a polytope, finds a set of occluders that is intersected by the set of lines that the polytope represents.

        The occluder finding is done using a ray-tracing operation, the ray(s) direction used to perform the sampling beeing either a "representative line" of the polytope,
        or the extremal stabbing line of the polytope.
        */
        bool collectAllOccluders(PluckerPolytope<P>* polytope, PluckerPolyhedron<P>* polyhedron, std::vector<VisibilitySilhouette*>& occluders, std::vector<P>& polytopeLines);

        /**@brief Attach a debugger to the query for visual inspection */
        void attachVisualisationDebugger(VisualDebugger* aDebugger);

        /**@brief Compute a visibility query to determine if the two polygons given as input are mutually visible or not */
        VisibilityResult arePolygonsVisible(const float* vertices0, size_t numVertices0, const float* vertices1, size_t numVertices1);
        /*
        const std::unordered_set<PluckerPolytope<P>*>& getPolytopes(VisibilitySilhouette* silhouette)
        {
            return mSilhouetteToPolytopeDictionary.find(silhouette)->second;
        }

        bool hasSilhouette(PluckerPolytope<P>* polytope)
        {
            return mPolytopeToSilhouetteDictionary.find(polytope) != mPolytopeToSilhouetteDictionary.end();
        }
        */
       
        
        /**@brief Find the next edge to be processed by the query
        */
        bool findNextEdge(size_t& aSilhouetteEdgeIndex, VisibilitySilhouette * &aSilhouette, PluckerPolytope<P> * polytope, const std::string & occlusionTreeNodeSymbol);


        void setApproximateNormal(const MathVector3d & a)
        {
            mApproximateNormal = a;
        }

        const MathVector3d& getApproximateNormal() const
        {
            return mApproximateNormal;
        }
        bool isOccluded(PluckerPolytope<P> * polytope, PluckerPolyhedron<P>* polyhedron, const std::vector <VisibilitySilhouette*> & silhouettes, const std::vector<P> & polytopeLines);

        /**@brief Return the statistic collector */
        HelperStatisticCollector* getStatistic()
        {
            return &mStatistic;
        }

        /**@brief Output the computation statistic in the console */
        void displayStatistic()
        {
            getStatistic()->display();
        }
    private:

        /**@brief Create the initial source polygons from the polygons provided as input

        If the suport plane of one of the input polygon intersect the other polygon, we clip the intersected polygon using the equation of the support plane.
        After both input polygons have been clipped, we are sure that the supporting planes of the two source polygons do not intersect inside one of the source polygons
        */
        bool createInitialPolygons(const float* vertices0, size_t numVertices0, const float* vertices1, size_t numVertices1, bool normalization);

        QueryConfiguration mConfiguration;                     /**< @brief The configuration parameters of the query*/
        PluckerPolytopeComplex<P>* mComplex;                   /**< @brief The polytope complex encoding the occlusion tree*/
        GeometryOccluderSet* mScene;                                         /**< @brief The scene containing the triangle mesh occluders*/
        GeometryConvexPolygon* mQueryPolygon[2];               /**< @brief The two convex polygonal sources*/
        VisualDebugger* mDebugger;                             /**< @brief The visual debugging information of the query*/
        HelperStatisticCollector mStatistic;                   /**< @brief The statistics of the query*/
        SilhouetteProcessor* mSilhouetteProcessor;   /**< @brief The silhouette processor for silhouette optimization heuristic*/

        S mTolerance;                                           /**< @brief The tolerance parameter used during computation*/
        MathVector3d mApproximateNormal;
        /** @brief The links between the polytopes and the silhouettes*/
 //      std::unordered_map<VisibilitySilhouette*, std::unordered_set<PluckerPolytope<P>*>> mSilhouetteToPolytopeDictionary;

        SilhouetteContainer* mSilhouetteContainer;
        /** @brief The links between the silhouettes and the polytopes*/
    //    std::unordered_map<PluckerPolytope<P>*, std::unordered_set<VisibilitySilhouette*>> mPolytopeToSilhouetteDictionary;
    };

    template<class P, class S>
    VisibilityExactQuery<P, S>::VisibilityExactQuery(GeometryOccluderSet * aScene, const QueryConfiguration & aConfiguration, S aTolerance)
        : mScene(aScene),
        mConfiguration(aConfiguration),
        mDebugger(nullptr)
    {
        mComplex = new PluckerPolytopeComplex<P>();

        {       
            HelperScopedTimer timer(getStatistic(), SILHOUETTE_PROCESSING);
            mSilhouetteProcessor = new SilhouetteProcessor(&mStatistic);
        }
        mQueryPolygon[0] = nullptr;
        mQueryPolygon[1] = nullptr;

        mTolerance = aTolerance;
 #if EMBREE       
        if (aConfiguration.useEmbree)
        {
            mSilhouetteContainer = new SilhouetteContainerEmbree();
        }
        else
#endif
        {
            mSilhouetteContainer = new SilhouetteContainer();
        }
    }

    template<class P, class S>
    VisibilityExactQuery<P, S>::~VisibilityExactQuery()
    {
        delete mComplex;
        delete mSilhouetteProcessor;
        delete mQueryPolygon[0];
        delete mQueryPolygon[1];
        delete mSilhouetteContainer;
     }

    template<class P, class S>
    bool VisibilityExactQuery<P, S>::createInitialPolygons(const float* vertices0, size_t numVertices0, const float* vertices1, size_t numVertices1, bool)
    {
        GeometryConvexPolygon myQuery0(vertices0, numVertices0);
        GeometryConvexPolygon myQuery1(vertices1, numVertices1);

        MathVector3d g1 = MathGeometry::getGravityCenter(myQuery0);
        MathVector3d g2 = MathGeometry::getGravityCenter(myQuery1);

        MathVector3d myN = g2 - g1;

        if (myQuery0.getVertexCount() < 3 || myQuery1.getVertexCount() < 3)
        {

            if (myQuery0.getVertexCount() == 2)
            {
                myQuery0.setPlane(MathGeometry::computePlaneFromApproximateNormal(myQuery0, myN));
            }
            else if (myQuery0.getVertexCount() == 1)
            {
                MathVector3d myNormal = myN;
                myNormal.normalize();
                myQuery0.setPlane(MathPlane3d(myNormal.x, myNormal.y, myNormal.z, -myNormal.dot(myQuery0.getVertex(0))));
            }

            if (myQuery1.getVertexCount() == 2)
            {
                myN.x = -myN.x; myN.y = -myN.y; myN.z = -myN.z;
                myQuery1.setPlane(MathGeometry::computePlaneFromApproximateNormal(myQuery1, myN));
            }
            else if (myQuery1.getVertexCount() == 1)
            {
                myN.x = -myN.x; myN.y = -myN.y; myN.z = -myN.z;
                MathVector3d myNormal = myN;
                myNormal.normalize();
                myQuery1.setPlane(MathPlane3d(myNormal.x, myNormal.y, myNormal.z, -myNormal.dot(myQuery1.getVertex(0))));
            }
        }
        setApproximateNormal(myN);

        if (!MathGeometry::clipWithGardBand(myQuery0, myQuery1.getPlane(), MathArithmetic<double>::GuardBandClipping()))
        {
            std::cout << "Clipping failure" << std::endl;
            return false;
        }
        V_ASSERT(myQuery0.isValid());

        if (!MathGeometry::clipWithGardBand(myQuery1, myQuery0.getPlane(), MathArithmetic<double>::GuardBandClipping()))
        {
            std::cout << "Clipping failure" << std::endl;
            return false;
        }
        V_ASSERT(myQuery1.isValid());

        V_ASSERT(MathPredicates::getRelativePosition(myQuery0.getVertices(), myQuery1.getPlane()) == ON_POSITIVE_SIDE);
        V_ASSERT(MathPredicates::getRelativePosition(myQuery1.getVertices(), myQuery0.getPlane()) == ON_POSITIVE_SIDE);

        mQueryPolygon[0] = new GeometryConvexPolygon(myQuery0);
        mQueryPolygon[1] = new GeometryConvexPolygon(myQuery1);

        return true;
    }

    template<class P, class S>
    VisibilityResult VisibilityExactQuery<P, S>::arePolygonsVisible(const float* vertices0, size_t numVertices0, const float* vertices1, size_t numVertices1)
    {
        VisibilityResult result;

        HelperScopedTimer timer(&mStatistic, VISIBILITY_QUERY);

        if (mDebugger != nullptr)
        {
            mDebugger->get()->clear();
            mSilhouetteProcessor->attachVisualisationDebugger(mDebugger);
        }

        {
            HelperScopedTimer timerBuild(&mStatistic, POLYTOPE_BUILD);
            if (!createInitialPolygons(vertices0, numVertices0, vertices1, numVertices1, mConfiguration.hyperSphereNormalization))
            {
                std::cout << "Polygon initialization error" << std::endl;
                return FAILURE;
            }
        }

        if (mQueryPolygon[0]->getVertexCount() == 1 && mQueryPolygon[1]->getVertexCount() == 1)
        {
            std::set<SilhouetteMeshFace*> intersectedFaces;
            result = findSceneIntersection(mQueryPolygon[0]->getVertex(0), mQueryPolygon[1]->getVertex(0), intersectedFaces) ? HIDDEN : VISIBLE;
        }
        else
        {   
            {
                HelperScopedTimer timer(getStatistic(), SILHOUETTE_PROCESSING);
                //mScene->get()->restoreFacesGeometry(mScene);

                mSilhouetteProcessor->init(*mQueryPolygon[0], *mQueryPolygon[1]);
                extractAllSilhouettes();
            }
            {
                HelperScopedTimer timer(getStatistic(), RAY_INTERSECTION);
                mSilhouetteContainer->prepare();
            }

            {
                HelperScopedTimer timerBuild(&mStatistic, POLYTOPE_BUILD);

                PluckerPolytopeBuilder<P, S> builder(mConfiguration.hyperSphereNormalization, mTolerance);
                PluckerPolytope<P>* myPolytope = builder.build(*mQueryPolygon[0], *mQueryPolygon[1], getComplex()->getPolyhedron());
                getComplex()->setRoot(myPolytope);
            }
            VisibilitySolver<P, S>* solver;

            solver = new VisibilityApertureFinder<P, S>(this, mConfiguration.hyperSphereNormalization, mTolerance, mConfiguration.performFastTest, mConfiguration.detectApertureOnly);

            if (mDebugger != nullptr)
            {
                solver->attachVisualisationDebugger(mDebugger);
            }
            result = solver->resolve();

            delete solver;
        }

        return result;
    }

    template<class P, class S>
    void VisibilityExactQuery<P, S>::attachVisualisationDebugger(VisualDebugger * aDebugger)
    {
        mDebugger = aDebugger;
    }

    template<class P, class S>
    bool VisibilityExactQuery<P, S>::findNextEdge(size_t& aSilhouetteEdgeIndex, VisibilitySilhouette * &aSilhouette, PluckerPolytope<P> * aPolytope, const std::string & occlusionTreeNodeSymbol)
    {
        HelperScopedTimer timer(getStatistic(), OCCLUDER_TREATMENT);
        std::ofstream& debugOutput = mDebugger->get()->getDebugOutput();
        bool hasDebugOutput = mDebugger->get()->hasDebugOutput();

        V_LOG(debugOutput, "VisibilityExactQuery<P, S>::findTheBestValidEdge BEGIN", occlusionTreeNodeSymbol);

        const std::unordered_set<VisibilitySilhouette*>& mySilhouettes = mSilhouetteContainer->getSilhouettes();

        double myScore = 1e32;
        bool found = false;
        
        VisibilitySilhouette* mySilhouette = nullptr;

        MathPlane3d aPlane0 = getQueryPolygon(0)->getPlane();
   
        const MathPlane3d& myPlane = getQueryPolygon(0)->getPlane();
        for (auto iter = mySilhouettes.begin(); iter != mySilhouettes.end(); iter++)
        {
            VisibilitySilhouette* s = (*iter);

            auto& edges = s->getEdges();

            for (size_t silhouetteEdgeIndex = 0; silhouetteEdgeIndex < edges.size(); silhouetteEdgeIndex++)
            {
                SilhouetteEdge& edge = edges[silhouetteEdgeIndex];
                if (edge.mIsActive)
                {
                    /*
                    if (edge.mScore < 0.0)
                    {
                        ConnectedMeshElement* f = edge.mFace;
                        size_t edgeIndex = edge.mEdgeIndex;
                        MathVector2i ei = f->getEdge(edgeIndex);
                        edge.mScore = aPlane0.dot(convert<MathVector3d>(f->getVertex(ei.x)));
                    }
                    if (edge.mScore < myScore)
                    */
                    {
                        found = true;
                        myScore = edge.mScore;
                        aSilhouetteEdgeIndex = silhouetteEdgeIndex;
                        mySilhouette = s;
                    }
                }
                if (found)
                    break;
            }

            if (found)
                break;
        }

        if (found)
        {
            std::stringstream ss;
            if (hasDebugOutput)
                ss << "VisibilityExactQuery<P, S>::findTheBestValidEdge END return True (Edge found)";
            V_LOG(debugOutput, ss.str(), occlusionTreeNodeSymbol);
            aSilhouette = mySilhouette;

        }
        else
        {
            V_LOG(debugOutput, "VisibilityExactQuery<P, S>::findTheBestValidEdge END return False (No edge found)", occlusionTreeNodeSymbol);
        }
        return found;
    }

    template<class P, class S>
    bool VisibilityExactQuery<P, S>::isOccluded(PluckerPolytope<P>* polytope, PluckerPolyhedron<P>* polyhedron, const std::vector <VisibilitySilhouette*> & aSilhouettes, const std::vector<P> & polytopeLines)
    {   
        return SilhouetteContainer::isOccluded(polytope, polyhedron, aSilhouettes,polytopeLines,mTolerance);
    }

    template<class P, class S>
    bool VisibilityExactQuery<P, S>::findSceneIntersection(const MathVector3d & aBegin, const MathVector3d & anEnd, std::set<SilhouetteMeshFace*> & anIntersectedFaces, PluckerPolytope<P> * aPolytope)
    {
        MathVector3d myBegin = aBegin;
        MathVector3d myDir = anEnd - aBegin;
        double myMax = myDir.normalize();

        Ray myRay;
        myRay.org[0] = (float)myBegin.x;
        myRay.org[1] = (float)myBegin.y;
        myRay.org[2] = (float)myBegin.z;
        myRay.dir[0] = (float)myDir.x;
        myRay.dir[1] = (float)myDir.y;
        myRay.dir[2] = (float)myDir.z;
        myRay.tfar = (float)myMax;
        myRay.tnear = 0;

        if (mDebugger != nullptr)
        {
            mDebugger->get()->addSamplingLine(convert<MathVector3f>(aBegin), convert<MathVector3f>(anEnd));
        }

        bool intersect = false;

        {
            HelperScopedTimer timer(getStatistic(), RAY_INTERSECTION);
            getStatistic()->inc(RAY_COUNT);
            
            intersect = mSilhouetteContainer->intersect(&myRay);
        }

        if (intersect)
        {
            for (size_t i = 0; i < myRay.mPrimitiveIds.size(); i++)
            {
                //:TODO: replace wihth accessors
                size_t myPrimitive = myRay.mPrimitiveIds[i];
                size_t geometryId = myRay.mGeometryIds[i];

                std::vector<SilhouetteMeshFace>* myFaces = mScene->getOccluderConnectedFaces(geometryId);
                SilhouetteMeshFace* myFace = &(*myFaces)[myPrimitive];
                anIntersectedFaces.insert(myFace);
            }
            return true;
        }

        if (mDebugger != nullptr)
        {
            mDebugger->get()->addStabbingLine(convert<MathVector3f>(aBegin), convert<MathVector3f>(anEnd));
        }

        return false;
    }

    template<class P, class S>
    void VisibilityExactQuery<P, S>::extractAllSilhouettes()
    {
        for (size_t geometryId = 0; geometryId < mScene->getOccluderCount(); geometryId++)
        {
            std::vector<VisibilitySilhouette*> silhouettes;
            std::vector<SilhouetteMeshFace>* myFaces = mScene->getOccluderConnectedFaces(geometryId);

            mSilhouetteProcessor->extractSilhouette(geometryId, *myFaces, mConfiguration.silhouetteOptimization, silhouettes);
            
            for (auto s:silhouettes)
            {
                mSilhouetteContainer->addSilhouette(s);
            }
            
        }
    }

    template<class P, class S>
    bool VisibilityExactQuery<P, S>::collectAllOccluders(PluckerPolytope<P> * aPolytope, PluckerPolyhedron<P> * polyhedron, std::vector<VisibilitySilhouette*> & occluders, std::vector<P> & polytopeLines)
    {
        {   
            if (mConfiguration.representativeLineSampling)
            {
                HelperScopedTimer timer(getStatistic(), STABBING_LINE_EXTRACTION);

                P myRepresentativeLine = MathGeometry::computeRepresentativeLine<P>(aPolytope, polyhedron, mTolerance);
                if (mConfiguration.hyperSphereNormalization)
                    myRepresentativeLine = myRepresentativeLine.getNormalized();
                polytopeLines.push_back(myRepresentativeLine);

                aPolytope->setRepresentativeLine(myRepresentativeLine);
            }
            else
            {
                HelperScopedTimer timer(getStatistic(), STABBING_LINE_EXTRACTION);

                if (aPolytope->getExtremalStabbingLinesCount() == 0)
                {
                    aPolytope->computeExtremalStabbingLines(polyhedron, mTolerance);
                }
                for (size_t i = 0; i < aPolytope->getExtremalStabbingLinesCount(); i++)
                {
                    polytopeLines.push_back(aPolytope->getExtremalStabbingLine(i));
                }
            }
        }

        const MathPlane3d& aPlane0 = getQueryPolygon(0)->getPlane();
        const MathPlane3d& aPlane1 = getQueryPolygon(1)->getPlane();

        for (size_t i = 0; i < polytopeLines.size(); i++)
        {
            std::pair<MathVector3d, MathVector3d> line = MathGeometry::getBackTo3D(polytopeLines[i], aPlane0, aPlane1);
            // V_ASSERT(MathGeometry::isPointInsidePolygon(*getQueryPolygon(0), lines[i].first, MathArithmetic<double>::Tolerance()));
            // V_ASSERT(MathGeometry::isPointInsidePolygon(*getQueryPolygon(1), lines[i].second, MathArithmetic<double>::Tolerance()));
            std::set<SilhouetteMeshFace*> intersectedFaces;
            bool hit = findSceneIntersection(line.first, line.second, intersectedFaces, aPolytope);
       //     if (mConfiguration.useEmbree )
         //        hit = hit || findSceneIntersection(line.second, line.first, intersectedFaces, aPolytope);
         
            if (!hit)
            {
                return false;
            }

            for (auto myFace : intersectedFaces)
            {
                VisibilitySilhouette* s = mSilhouetteProcessor->findSilhouette(myFace);
                //   V_ASSERT(s);
                if (s)
                    occluders.push_back(s);
            }
        }
        return true;
    }
}
