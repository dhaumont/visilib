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

#include <unordered_map>
#include <unordered_set>
#include <random>
#include <cmath>
#include "geometry_convex_polygon.h"
#include "silhouette_mesh_face.h"
#include "geometry_occluder_set.h"
#include "helper_statistic_collector.h"
#include "math_geometry.h"
#include "math_predicates.h"
#include "math_plucker_2.h"
#include "math_vector_3.h"
#include "plucker_polytope.h"
#include "plucker_polytope_builder.h"
#include "plucker_polytope_complex.h"
#include "visibility_aperture_finder.h"
#include "visibility_aggressive_solver.h"
#include "visibility_sequential_solver.h"
#include "silhouette_container.h"
#include "silhouette_processor.h"
#include "visilib.h"
#include "visilib_core.h"
#include "silhouette_container_embree.h"
#include "shadow_discontinuity_extractor.h"
#include "visibility_event.h"
#include "visibility_event_registry.h"

namespace visilib
{
    class HelperDebugVisualisation;
    class HelperStatisticCollector;

    class GeometryConvexPolygon;
    class SilhouetteMeshFace;
    class GeometryOccluderSet;
    struct VisibilityRay;
    class Silhouette;
    class SilhouetteProcessor;

    template <class S>
    class PluckerPolytope;

    /**@brief Exact visibility query interface class

    This interface is the protoype definition of an exact visibility query between two source polygons
    */
    class VisibilityExactQuery
    {
    public:
        virtual ~VisibilityExactQuery() {};

        virtual void attachVisualisationDebugger(HelperVisualDebugger* aDebugger) = 0;
        virtual VisibilityResult arePolygonsVisible(const float* vertices0, size_t numVertices0, const float* vertices1, size_t numVertices1) = 0;
        virtual void displayStatistic() = 0;
        virtual void getTimingStats(std::map<std::string, double> & stats) = 0;
        virtual void getCounts(std::map<std::string, int> & stats) = 0;

        /**@brief Return the measurements of the query, the values that are neither tallies nor times

        Not a pure virtual: a query that measures nothing has nothing to say about it, and an
        implementation written before the question existed is not made wrong by it being asked. */
        virtual void getMeasures(std::map<std::string, double> & /*measures*/) {}

        /**@brief Return the shadow discontinuities traced on the receiver (query polygon 1) by the light (query polygon 0)
        during the last query. Empty unless VisibilityExactQueryConfiguration::computeShadowDiscontinuities is set.*/
        const std::vector<ShadowDiscontinuityArc>& getShadowDiscontinuities() const
        {
            return mShadowDiscontinuities;
        }

    protected:
        std::vector<ShadowDiscontinuityArc> mShadowDiscontinuities;   /**< @brief The discontinuity arcs of the last query*/
    };


    /**@brief Exact Visibility solver implementation class

    This class manages the computation of exact visibility query between two polygons through a set of occluders contained in a scene.
    */
    template<class P, class S>
    class VisibilityExactQuery_ : public VisibilityExactQuery
    {
    public:

        VisibilityExactQuery_(GeometryOccluderSet* aScene, const VisibilityExactQueryConfiguration& aConfiguration, S tolerance);

        virtual ~VisibilityExactQuery_();

        /**@brief Whether the query gathers the shadow discontinuities of its visible cells*/
        bool isComputingShadowDiscontinuities() const
        {
            return mShadowExtractor != nullptr;
        }

        /**@brief Record that the hyperplane anIndex of the polyhedron carries the occluder edge [a, b], aFlat when it joins two coplanar faces*/
        void registerOccluderEdge(size_t anIndex, const MathVector3d& a, const MathVector3d& b, bool aFlat)
        {
            if (mShadowExtractor == nullptr)
            {
                return;
            }
            VisibilityFeature feature;
            feature.mKind = VisibilityFeature::OCCLUDER_EDGE;
            feature.mEndpoints[0] = a;
            feature.mEndpoints[1] = b;
            feature.mRealEndpoint[0] = true;
            feature.mRealEndpoint[1] = true;
            feature.mSynthetic = false;
            feature.mFlat = aFlat;
            mEventRegistry.registerHyperplane(anIndex, feature);
        }

        /**@brief Gather the shadow discontinuities of a cell certified visible. The cell must be alive.*/
        void addVisibleCell(PluckerPolyhedron<P>* aPolyhedron, PluckerPolytope<P>* aPolytope)
        {
            if (mShadowExtractor != nullptr)
            {
                mShadowExtractor->addVisibleCell(aPolyhedron, aPolytope, mTolerance);
            }
        }

        /**@brief Return the source polygons at the given index */
        GeometryConvexPolygon* getQueryPolygon(size_t i)
        {
            return mQueryPolygon[i];
        }

        /**@brief Return the occluder scene */
        GeometryOccluderSet* getScene()
        {
            return mScene;
        }

        /**@brief Return the polytope complex encoding the occlusion tree */
        PluckerPolytopeComplex<P>* getComplex()
        {
            return mComplex;
        }

        /**@brief Performs the intersection between a segment and the geometry of the scene.

        When an intersection is found, the silhouette are extracted from the intersected face and the silhouette is linked to the polytope
        The links is used for occluder selection.

        With a radius at each end of the segment the test is not the thin line but the beam
        capsule that spans it — every line joining a point of the disk at one end to a point
        of the disk at the other. A cast of that body finds every face any line of the beam
        can touch: an occluder selection made with it cannot miss one the segment itself
        grazes past.
        */
        bool findSceneIntersection(const MathVector3d& aBegin, const MathVector3d& anEnd, std::set<SilhouetteMeshFace*>& intersectedFaces, double radiusAtStart = 0, double radiusAtEnd = 0);
        
        bool hasSceneIntersection(const MathVector3d& aBegin, const MathVector3d& anEnd);

        void extractAllSilhouettes();

        /**@brief Given a polytope, finds a set of occluders that is intersected by the set of lines that the polytope represents.

        The occluder finding is done using a ray-tracing operation, the ray(s) direction used to perform the sampling beeing either a "representative line" of the polytope,
        or the extremal stabbing line of the polytope.
        */
        int collectAllOccluders(PluckerPolytope<P>* polytope, PluckerPolyhedron<P>* polyhedron, std::vector<Silhouette*>& occluders, std::vector<P>& polytopeLines);

        /**@brief Compute the beam capsule that encloses every line of a polytope cell.

        The extremal stabbing lines of the cell pierce the two query planes at the points
        that bound the piercing region of the cell: the enclosing disk of each end's points
        gives one ball of the capsule, and every line of the cell — a segment joining a
        point of one region to a point of the other — runs through the body between them.

        The two centers are the mean of each end's points and the two radii the greatest
        distance to them, the representative line's own piercing points counted among them;
        a smallest enclosing circle would tighten each end further and is left as a later
        refinement.

        A stabbing line piercing far outside its polygon is an artifact of the arithmetic —
        a line near-parallel to the plane has its piercing point anywhere — and is left out
        of both disks, while every real vertex of the cell is checked against them and the
        radius grown to hold one that escapes: the beam's purpose is to miss no occluder,
        so it errs by width and never by measure.

        Returns false when the cell has nothing to bound the region with; the caller then
        abstains, as it did before any beam was asked.
        */
        bool computeOccluderSearchBeam(PluckerPolytope<P>* polytope, PluckerPolyhedron<P>* polyhedron, const P& representativeLine, MathVector3d& beamStart, MathVector3d& beamEnd, double& radiusStart, double& radiusEnd);

        /**@brief Check that the representative line of a leaf polytope is a valid witness of visibility: the line must intersect the two source polygons and must not be blocked by the scene.

        Returns 1 if the representative line is a valid, unoccluded witness, 0 if it is occluded, and -1 if the sample is invalid (its endpoints are not inside both source polygons).
        */
        int validateLeafPolytope(PluckerPolytope<P>* polytope, PluckerPolyhedron<P>* polyhedron);

        /**@brief Attach a debugger to the query for visual inspection */
        void attachVisualisationDebugger(HelperVisualDebugger* aDebugger);

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
        bool findNextEdge(size_t& aSilhouetteEdgeIndex, Silhouette * &aSilhouette, PluckerPolytope<P> * polytope, const std::string & occlusionTreeNodeSymbol);


        const std::vector<Silhouette*>& getSilhouettes()
        {
            return mSilhouetteContainer->getSilhouettes();
        }

        void setApproximateNormal(const MathVector3d & a)
        {
            mApproximateNormal = a;
        }

        const MathVector3d& getApproximateNormal() const
        {
            return mApproximateNormal;
        }
        bool isOccluded(PluckerPolytope<P> * polytope, PluckerPolyhedron<P>* polyhedron, const std::vector <Silhouette*> & silhouettes, const std::vector<P> & polytopeLines);

        /**@brief Return the configuration the query is running with */
        const VisibilityExactQueryConfiguration& getConfiguration() const
        {
            return mConfiguration;
        }

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

        void getTimingStats(std::map<std::string, double> & stats)
        {
            return mStatistic.getTimingStats(stats);
        }

        void getCounts(std::map<std::string, int> & stats)
        {
            return mStatistic.getCounts(stats);
        }

        void getMeasures(std::map<std::string, double> & stats)
        {
            return mStatistic.getMeasures(stats);
        }
    private:

        /**@brief Create the initial source polygons from the polygons provided as input

        If the suport plane of one of the input polygon intersect the other polygon, we clip the intersected polygon using the equation of the support plane.
        After both input polygons have been clipped, we are sure that the supporting planes of the two source polygons do not intersect inside one of the source polygons
        */
        bool createInitialPolygons(const float* vertices0, size_t numVertices0, const float* vertices1, size_t numVertices1, bool normalization);

        VisibilityExactQueryConfiguration mConfiguration;                     /**< @brief The configuration parameters of the query*/
        PluckerPolytopeComplex<P>* mComplex;                   /**< @brief The polytope complex encoding the occlusion tree*/
        GeometryOccluderSet* mScene;                                         /**< @brief The scene containing the triangle mesh occluders*/
        GeometryConvexPolygon* mQueryPolygon[2];               /**< @brief The two convex polygonal sources*/
        HelperVisualDebugger* mDebugger;                             /**< @brief The visual debugging information of the query*/
        HelperStatisticCollector mStatistic;                   /**< @brief The statistics of the query*/
        SilhouetteProcessor* mSilhouetteProcessor;   /**< @brief The silhouette processor for silhouette optimization heuristic*/

        S mTolerance;                                           /**< @brief The tolerance parameter used during computation*/
        MathVector3d mApproximateNormal;
        /** @brief The links between the polytopes and the silhouettes*/
 //      std::unordered_map<VisibilitySilhouette*, std::unordered_set<PluckerPolytope<P>*>> mSilhouetteToPolytopeDictionary;

        SilhouetteContainer* mSilhouetteContainer;

        VisibilityEventRegistry mEventRegistry;                          /**< @brief The edges the hyperplanes of the polyhedron carry*/
        ShadowDiscontinuityExtractor<P, S>* mShadowExtractor;            /**< @brief The gathering of the shadow discontinuities, during a query that computes them*/
        /** @brief The links between the silhouettes and the polytopes*/
    //    std::unordered_map<PluckerPolytope<P>*, std::unordered_set<VisibilitySilhouette*>> mPolytopeToSilhouetteDictionary;
    };

    template<class P, class S>
    VisibilityExactQuery_<P, S>::VisibilityExactQuery_(GeometryOccluderSet * aScene, const VisibilityExactQueryConfiguration & aConfiguration, S aTolerance)
        : mScene(aScene),
        mConfiguration(aConfiguration),
        mDebugger(nullptr),
        mShadowExtractor(nullptr)
    {
        mComplex = new PluckerPolytopeComplex<P>();

        mStatistic.setDisplaySplitOccurence(aConfiguration.displaySplitOccurence);

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
    VisibilityExactQuery_<P, S>::~VisibilityExactQuery_()
    {
        delete mComplex;
        delete mSilhouetteProcessor;
        delete mQueryPolygon[0];
        delete mQueryPolygon[1];
        delete mSilhouetteContainer;
        delete mShadowExtractor;
     }

    template<class P, class S>
    bool VisibilityExactQuery_<P, S>::createInitialPolygons(const float* vertices0, size_t numVertices0, const float* vertices1, size_t numVertices1, bool)
    {
        GeometryConvexPolygon myQuery0(vertices0, numVertices0);
        GeometryConvexPolygon myQuery1(vertices1, numVertices1);

        MathVector3d g1 = MathGeometry::getGravityCenter(myQuery0);
        MathVector3d g2 = MathGeometry::getGravityCenter(myQuery1);

        MathVector3d myN = g2 - g1;
        myN.normalize();
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
    VisibilityResult VisibilityExactQuery_<P, S>::arePolygonsVisible(const float* vertices0, size_t numVertices0, const float* vertices1, size_t numVertices1)
    {
        VisibilityResult result = UNKNOWN;

        HelperScopedTimer timer(&mStatistic, VISIBILITY_QUERY);

        if (mDebugger != nullptr)
        {
            mDebugger->clear();
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

        mShadowDiscontinuities.clear();
        mEventRegistry.clear();
        delete mShadowExtractor;
        mShadowExtractor = nullptr;
        if (mConfiguration.computeShadowDiscontinuities)
        {
            /* The shadows are traced on the plane of the receiver: a point or a segment has none. The discontinuities are
            read from every visible cell, which only the aperture finder in its complete mode produces. */
            if (mQueryPolygon[1]->getVertexCount() < 3)
            {
                std::cerr << "Shadow discontinuities: the receiver (query polygon 1) must be a polygon" << std::endl;
            }
            else if (mConfiguration.solverType != VisibilityExactQueryConfiguration::EXACT_APERTURE_FINDER || mConfiguration.detectApertureOnly)
            {
                std::cerr << "Shadow discontinuities: the aperture finder with detectApertureOnly off is required" << std::endl;
            }
            else
            {
                mShadowExtractor = new ShadowDiscontinuityExtractor<P, S>(&mEventRegistry, *mQueryPolygon[0], *mQueryPolygon[1]);
            }
        }

        if (mQueryPolygon[0]->getVertexCount() == 1 && mQueryPolygon[1]->getVertexCount() == 1)
        {
            /* The rays are cast against the silhouette container: it has to be filled first, or the segment
            meets nothing and every pair of points is answered visible. */
            {
                HelperScopedTimer timer(getStatistic(), SILHOUETTE_PROCESSING);
                mSilhouetteProcessor->init(*mQueryPolygon[0], *mQueryPolygon[1]);
                extractAllSilhouettes();
            }
            {
                HelperScopedTimer timer(getStatistic(), RAY_INTERSECTION);
                mSilhouetteContainer->prepare();
            }
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

                PluckerPolytopeBuilder<P, S> builder(mConfiguration.hyperSphereNormalization, mTolerance, mShadowExtractor != nullptr ? &mEventRegistry : nullptr);
                PluckerPolytope<P>* myPolytope = builder.build(*mQueryPolygon[0], *mQueryPolygon[1], getComplex()->getPolyhedron());
                getComplex()->appendPolytope(myPolytope);
            }
            VisibilitySolver<P, S>* solver;

            switch (mConfiguration.solverType)
            {
                case VisibilityExactQueryConfiguration::EXACT_APERTURE_FINDER:
                    solver = new VisibilityApertureFinder<P, S>(this, mConfiguration.hyperSphereNormalization, mTolerance, mConfiguration.detectApertureOnly);
                break;

                case VisibilityExactQueryConfiguration::MONTE_CARLO:
                    solver = new VisibilityAggressiveSolver<P, S>(this, mTolerance,
                                       mConfiguration.detectApertureOnly, 
                                       mConfiguration.minimumNormalizedApertureSize,
                                       mConfiguration.confidenceValue);
                break;                
                case VisibilityExactQueryConfiguration::EXACT_SEQUENTIAL_SOLVER:
                    solver = new VisibilitySequentialSolver<P, S>(this, mConfiguration.hyperSphereNormalization, mTolerance, mConfiguration.detectApertureOnly);
                break;
            }
            if (mDebugger != nullptr)
            {
                solver->attachVisualisationDebugger(mDebugger);
            }
            result = solver->resolve();
          
            delete solver;
        }

        if (mShadowExtractor != nullptr)
        {
            mShadowExtractor->getArcs(mShadowDiscontinuities);

            if (mDebugger != nullptr)
            {
                /* Each arc as a polyline sampled along its family, lifted off the receiver plane by a hair
                towards the light, so that it is drawn over the receiver rather than in it. */
                const double lift = 1e-3 * (mQueryPolygon[1]->getVertex(0) - MathGeometry::getGravityCenter(*mQueryPolygon[1])).getNorm();
                for (const ShadowDiscontinuityArc& arc : mShadowDiscontinuities)
                {
                    MathVector3d offset = arc.mFrame.getNormal() * lift;
                    if (arc.mFrame.getSignedDistance(MathGeometry::getGravityCenter(*mQueryPolygon[0])) < 0.0)
                    {
                        offset *= -1.0;
                    }
                    const std::vector<MathVector2d> polyline = arc.toPolyline(arc.mType == ShadowDiscontinuityArc::CONIC ? 32 : 2);
                    for (size_t i = 1; i < polyline.size(); i++)
                    {
                        mDebugger->addShadowDiscontinuity(convert<MathVector3f>(arc.toPoint3D(polyline[i - 1]) + offset),
                                                          convert<MathVector3f>(arc.toPoint3D(polyline[i]) + offset));
                    }
                }
            }
            delete mShadowExtractor;
            mShadowExtractor = nullptr;
        }

        return result;
    }

    template<class P, class S>
    void VisibilityExactQuery_<P, S>::attachVisualisationDebugger(HelperVisualDebugger * aDebugger)
    {
        mDebugger = aDebugger;
    }

    template<class P, class S>
    bool VisibilityExactQuery_<P, S>::findNextEdge(size_t& aSilhouetteEdgeIndex, Silhouette * &aSilhouette, PluckerPolytope<P> * aPolytope, const std::string & occlusionTreeNodeSymbol)
    {
        HelperScopedTimer timer(getStatistic(), OCCLUDER_TREATMENT);
#ifdef OUTPUT_DEBUG_FILE
        std::ofstream& debugOutput = mDebugger->getDebugOutput();
        V_LOG(debugOutput, "VisibilityExactQuery<P, S>::findTheBestValidEdge BEGIN", occlusionTreeNodeSymbol);
#endif

        const std::vector<Silhouette*>& mySilhouettes = mSilhouetteContainer->getSilhouettes();

        double myScore = 1e32;
        bool found = false;

        Silhouette* mySilhouette = nullptr;

        MathPlane3d aPlane0 = getQueryPolygon(0)->getPlane();

        const MathPlane3d& myPlane = getQueryPolygon(0)->getPlane();
        for (auto iter = mySilhouettes.begin(); iter != mySilhouettes.end(); iter++)
        {
            Silhouette* s = (*iter);

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
#ifdef OUTPUT_DEBUG_FILE
            std::stringstream ss;
            ss << "VisibilityExactQuery<P, S>::findTheBestValidEdge END return True (Edge found)";
            V_LOG(debugOutput, ss.str(), occlusionTreeNodeSymbol);
#endif
            aSilhouette = mySilhouette;

        }
        else
        {
#ifdef OUTPUT_DEBUG_FILE
            V_LOG(debugOutput, "VisibilityExactQuery<P, S>::findTheBestValidEdge END return False (No edge found)", occlusionTreeNodeSymbol);
#endif
        }
        return found;
    }

    template<class P, class S>
    bool VisibilityExactQuery_<P, S>::isOccluded(PluckerPolytope<P>* polytope, PluckerPolyhedron<P>* polyhedron, const std::vector <Silhouette*> & aSilhouettes, const std::vector<P> & polytopeLines)
    {
        return SilhouetteContainer::isOccluded(polytope, polyhedron, aSilhouettes, polytopeLines, mTolerance, mConfiguration.computeShadowDiscontinuities);
    }

    template<class P, class S>
    bool VisibilityExactQuery_<P, S>::findSceneIntersection(const MathVector3d & aBegin, const MathVector3d & anEnd, std::set<SilhouetteMeshFace*> & anIntersectedFaces, double radiusAtStart, double radiusAtEnd)
    {
        const bool thin = (radiusAtStart == 0.0 && radiusAtEnd == 0.0);

        /* Early stop traces no fat ray: it stops at the first aperture it can prove, and the
        body a wide cast searches with is the whole cell — the very thing it declined to
        examine. A wide cast reaching here is a caller that asked the finder for work it is
        switched off to do, so the beam is refused rather than traced: the answer stays the
        thin one, which is all early stop ever meant to have. */
        if (!thin && mConfiguration.detectApertureOnly)
        {
            V_ASSERT(false);
            return false;
        }

        MathVector3d myBegin = aBegin;
        MathVector3d myDir = anEnd - aBegin;
        double myMax = myDir.normalize();

        VisibilityRay myRay;
        myRay.org[0] = (float)myBegin.x;
        myRay.org[1] = (float)myBegin.y;
        myRay.org[2] = (float)myBegin.z;
        myRay.dir[0] = (float)myDir.x;
        myRay.dir[1] = (float)myDir.y;
        myRay.dir[2] = (float)myDir.z;
        myRay.tfar = (float)myMax;
        myRay.tnear = 0;

        if (mDebugger != nullptr && thin)
        {
            mDebugger->addSamplingLine(convert<MathVector3f>(aBegin), convert<MathVector3f>(anEnd));
        }

        bool intersect = false;

        {
            HelperScopedTimer timer(getStatistic(), RAY_INTERSECTION);
            getStatistic()->inc(RAY_COUNT);

            intersect = mSilhouetteContainer->intersect(&myRay, radiusAtStart, radiusAtEnd);
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

        if (mDebugger != nullptr && thin)
        {
            mDebugger->addStabbingLine(convert<MathVector3f>(aBegin), convert<MathVector3f>(anEnd));
        }

        return false;
    }

    template<class P, class S>
    bool VisibilityExactQuery_<P, S>::hasSceneIntersection(const MathVector3d & aBegin, const MathVector3d & anEnd)
    {
        MathVector3d myBegin = aBegin;
        MathVector3d myDir = anEnd - aBegin;
        double myMax = myDir.normalize();

        VisibilityRay myRay;
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
            mDebugger->addSamplingLine(convert<MathVector3f>(aBegin), convert<MathVector3f>(anEnd));
        }

        bool intersect = false;

        {
            HelperScopedTimer timer(getStatistic(), RAY_INTERSECTION);
            getStatistic()->inc(RAY_COUNT);

            intersect = mSilhouetteContainer->intersect(&myRay);
        }

        return intersect;
    }
    template<class P, class S>
    void VisibilityExactQuery_<P, S>::extractAllSilhouettes()
    {
        for (size_t geometryId = 0; geometryId < mScene->getOccluderCount(); geometryId++)
        {
            std::vector<Silhouette*> silhouettes;
            std::vector<SilhouetteMeshFace>* myFaces = mScene->getOccluderConnectedFaces(geometryId);

            mSilhouetteProcessor->extractSilhouette(geometryId, *myFaces, mConfiguration.silhouetteOptimization, silhouettes);

            for (auto s:silhouettes)
            {
                mSilhouetteContainer->addSilhouette(s);
            }

        }
    }

    template<class P, class S>
    int VisibilityExactQuery_<P, S>::collectAllOccluders(PluckerPolytope<P> * aPolytope, PluckerPolyhedron<P> * polyhedron, std::vector<Silhouette*> & occluders, std::vector<P> & polytopeLines)
    {
        RepresentativeLine<P> myRepresentative;
        {
             HelperScopedTimer timer(getStatistic(), STABBING_LINE_EXTRACTION);

             myRepresentative = MathGeometry::representativeLine<P, S>(aPolytope, polyhedron, mTolerance, mConfiguration.optimisedRepresentativeLine, &mStatistic);
        }

        P myRepresentativeLine = myRepresentative.line;

        if (mConfiguration.hyperSphereNormalization)
            myRepresentativeLine = myRepresentativeLine.getNormalized();

        polytopeLines.push_back(myRepresentativeLine);
        aPolytope->setRepresentativeLine(myRepresentativeLine);
        aPolytope->setRepresentativeLineMargin(myRepresentative.margin);

        /* Which occluders a cell still has to be cut by is not to be read from a line lying outside of it: the cell
        is left alone, and the search comes back to it with a smaller one. */
        if (mConfiguration.optimisedRepresentativeLine && !myRepresentative.insideCell && !myRepresentative.boundary)
        {
            return -1;
        }

        const GeometryConvexPolygon& polygon0 = *getQueryPolygon(0);
        const GeometryConvexPolygon& polygon1 = *getQueryPolygon(1);
        const MathPlane3d& aPlane0 = polygon0.getPlane();
        const MathPlane3d& aPlane1 = polygon1.getPlane();

        std::pair<MathVector3d, MathVector3d> centerLine = MathGeometry::getBackTo3D(myRepresentativeLine, aPlane0, aPlane1);

        bool isInsidePolygons = MathGeometry::isPointInsidePolygon(polygon0, centerLine.first, MathArithmetic<double>::Tolerance())
            && MathGeometry::isPointInsidePolygon(polygon1, centerLine.second, MathArithmetic<double>::Tolerance());
    
        if (!isInsidePolygons)
        {          
            return -1;
        }

        std::set<SilhouetteMeshFace*> intersectedFaces;
        const bool hit = findSceneIntersection(centerLine.first, centerLine.second, intersectedFaces);

        /* The thin line went through and saw nothing: what it crossed is not nothing. Cast
        beside it the beam capsule that holds every line of this cell — the wall a grazing
        aperture hides is a face the beam reaches — and take what it finds for the cell's
        occluders. The verdict of the thin line stands, and nought is what comes back: the
        cell is visible, and the gathered silhouettes are what the subdivision still has to
        peel the shadow off. */
        if (!hit && !mConfiguration.detectApertureOnly)
        {
            MathVector3d beamStart;
            MathVector3d beamEnd;
            double radiusStart = 0.0;
            double radiusEnd = 0.0;

            if (computeOccluderSearchBeam(aPolytope, polyhedron, myRepresentativeLine, beamStart, beamEnd, radiusStart, radiusEnd))
            {
                getStatistic()->inc(FAT_RAY_CAST_COUNT);

                const size_t facesBefore = intersectedFaces.size();
                if (findSceneIntersection(beamStart, beamEnd, intersectedFaces, radiusStart, radiusEnd))
                {
                    getStatistic()->add(FAT_RAY_FACE_COUNT, intersectedFaces.size() - facesBefore);
                }
            }
            else
            {
                /* Nothing to bound the cell with, nothing to answer the cell from: it is left
                undecided, as the pass always did when its line could not be trusted. */
                return -1;
            }
        }

        for (auto myFace : intersectedFaces)
        {
            Silhouette* s = mSilhouetteProcessor->findSilhouette(myFace);
            //   V_ASSERT(s);
            if (s)
                occluders.push_back(s);
        }
        return hit ? 1 : 0;
    }

    template<class P, class S>
    bool VisibilityExactQuery_<P, S>::computeOccluderSearchBeam(PluckerPolytope<P> * aPolytope, PluckerPolyhedron<P> * polyhedron, const P & representativeLine, MathVector3d & beamStart, MathVector3d & beamEnd, double & radiusStart, double & radiusEnd)
    {
        {
            HelperScopedTimer timer(getStatistic(), STABBING_LINE_EXTRACTION);

            if (aPolytope->getExtremalStabbingLinesCount() == 0)
            {
                aPolytope->computeExtremalStabbingLines(polyhedron, mTolerance);
            }
        }

        if (aPolytope->getExtremalStabbingLinesCount() == 0)
        {
            return false;
        }

        const GeometryConvexPolygon& polygon0 = *getQueryPolygon(0);
        const GeometryConvexPolygon& polygon1 = *getQueryPolygon(1);
        const MathPlane3d& aPlane0 = polygon0.getPlane();
        const MathPlane3d& aPlane1 = polygon1.getPlane();

        /* The safety disk of a polygon: a line of the cell pierces its polygon, so a
        piercing point fallen well outside stands no chance of being an extreme of the cell
        at all — the back to three dimensions of a stabbing line run near-parallel to the
        plane lands anywhere, and a point like that must bound nothing. */
        auto safetyDisk = [](const GeometryConvexPolygon& polygon)
        {
            const MathVector3d center = MathGeometry::getGravityCenter(polygon);
            double radius = 0.0;
            for (size_t vertex = 0; vertex < polygon.getVertexCount(); vertex++)
            {
                MathVector3d toVertex = convert<MathVector3d>(polygon.getVertex(vertex)) - center;
                const double distance = toVertex.getNorm();
                if (distance > radius)
                {
                    radius = distance;
                }
            }
            return std::make_pair(center, radius * 1.25 + MathArithmetic<double>::Tolerance());
        };

        const std::pair<MathVector3d, double> disk0 = safetyDisk(polygon0);
        const std::pair<MathVector3d, double> disk1 = safetyDisk(polygon1);

        auto piercesWithin = [](const MathVector3d& point, const std::pair<MathVector3d, double>& disk)
        {
            if (!std::isfinite(point.x) || !std::isfinite(point.y) || !std::isfinite(point.z))
            {
                return false;
            }
            MathVector3d toPoint = point - disk.first;
            return toPoint.getSquaredNorm() <= disk.second * disk.second;
        };

        std::vector<MathVector3d> ends0;
        std::vector<MathVector3d> ends1;

        for (size_t i = 0; i < aPolytope->getExtremalStabbingLinesCount(); i++)
        {
            const std::pair<MathVector3d, MathVector3d> line = MathGeometry::getBackTo3D(aPolytope->getExtremalStabbingLine(i), aPlane0, aPlane1);

            if (piercesWithin(line.first, disk0) && piercesWithin(line.second, disk1))
            {
                ends0.push_back(line.first);
                ends1.push_back(line.second);
            }
            else
            {
                V_TRACE("extremal stabbing line pierces outside its polygon: left out of the beam");
            }
        }

        if (ends0.empty())
        {
            return false;
        }

        auto enclose = [](const std::vector<MathVector3d>& points, MathVector3d& center, double& radius)
        {
            center = MathVector3d();
            for (size_t i = 0; i < points.size(); i++)
            {
                center += points[i];
            }
            center *= 1.0 / (double)points.size();

            radius = 0.0;
            for (size_t i = 0; i < points.size(); i++)
            {
                MathVector3d toPoint = points[i] - center;
                const double distance = toPoint.getNorm();
                if (distance > radius)
                {
                    radius = distance;
                }
            }
        };

        /* The mean and the greatest distance per end: a smallest enclosing circle would
        tighten both ends further, and is the refinement left standing for later. */
        enclose(ends0, beamStart, radiusStart);
        enclose(ends1, beamEnd, radiusEnd);

        /* The representative line is a line of the cell: its piercing points, which the
        caller held inside the polygons, bound as much as any stabbing line does. */
        const std::pair<MathVector3d, MathVector3d> centerLine = MathGeometry::getBackTo3D(representativeLine, aPlane0, aPlane1);
        MathVector3d toStart = centerLine.first - beamStart;
        radiusStart = std::max(radiusStart, toStart.getNorm());
        MathVector3d toEnd = centerLine.second - beamEnd;
        radiusEnd = std::max(radiusEnd, toEnd.getNorm());

        /* The guard the beam keeps on itself: the disk of an end must hold the piercing
        points of every real line of the cell, not only of the stabbing lines believed to
        bound them. The cell's vertices that are real lines are such lines, and they are at
        hand: the radius that fails to hold one grows to hold it, and says so. */
        for (size_t vertexIndex : aPolytope->getVertices())
        {
            if (polyhedron->getQuadricRelativePosition(vertexIndex) != ON_BOUNDARY)
            {
                continue;
            }

            const std::pair<MathVector3d, MathVector3d> line = MathGeometry::getBackTo3D(polyhedron->get(vertexIndex), aPlane0, aPlane1);

            if (piercesWithin(line.first, disk0) && (line.first - beamStart).getNorm() > radiusStart)
            {
                V_TRACE("beam grown at its start for a cell vertex the stabbing lines left outside");
                radiusStart = (line.first - beamStart).getNorm();
            }
            if (piercesWithin(line.second, disk1) && (line.second - beamEnd).getNorm() > radiusEnd)
            {
                V_TRACE("beam grown at its end for a cell vertex the stabbing lines left outside");
                radiusEnd = (line.second - beamEnd).getNorm();
            }
        }

        /* A breath of room over the exact enclosing: the beam is cast in floats, and the
        face it must not miss is often the one it grazes. */
        radiusStart = radiusStart * 1.0001 + MathArithmetic<double>::Tolerance();
        radiusEnd = radiusEnd * 1.0001 + MathArithmetic<double>::Tolerance();

        return true;
    }

    template<class P, class S>
    int VisibilityExactQuery_<P, S>::validateLeafPolytope(PluckerPolytope<P> * aPolytope, PluckerPolyhedron<P> * polyhedron)
    {
        RepresentativeLine<P> myRepresentative;
        {
            HelperScopedTimer timer(getStatistic(), STABBING_LINE_EXTRACTION);

            myRepresentative = MathGeometry::representativeLine<P, S>(aPolytope, polyhedron, mTolerance, mConfiguration.optimisedRepresentativeLine, &mStatistic);
        }

        P myRepresentativeLine = myRepresentative.line;

        if (mConfiguration.hyperSphereNormalization)
            myRepresentativeLine = myRepresentativeLine.getNormalized();

        /* A leaf is certified visible from a line that is a real line of the cell: one outside of it says nothing,
        and saying so here keeps the answer short of a wrong one. */
        if (mConfiguration.optimisedRepresentativeLine && !myRepresentative.insideCell && !myRepresentative.boundary)
        {
            return -1;
        }

        const GeometryConvexPolygon& polygon0 = *getQueryPolygon(0);
        const GeometryConvexPolygon& polygon1 = *getQueryPolygon(1);
        const MathPlane3d& aPlane0 = polygon0.getPlane();
        const MathPlane3d& aPlane1 = polygon1.getPlane();

        std::pair<MathVector3d, MathVector3d> centerLine = MathGeometry::getBackTo3D(myRepresentativeLine, aPlane0, aPlane1);

        bool isInsidePolygons = MathGeometry::isPointInsidePolygon(polygon0, centerLine.first, MathArithmetic<double>::Tolerance())
            && MathGeometry::isPointInsidePolygon(polygon1, centerLine.second, MathArithmetic<double>::Tolerance());

        if (!isInsidePolygons)
        {
            return -1;
        }

        std::set<SilhouetteMeshFace*> intersectedFaces;
        bool hit = findSceneIntersection(centerLine.first, centerLine.second, intersectedFaces);

        return hit ? 0 : 1;
    }

}
