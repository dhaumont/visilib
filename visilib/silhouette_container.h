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
#include "visilib.h"
#include "math_vector_2.h"
#include "math_vector_3.h"
namespace visilib
{

 /** @brief
   Store a set of silhouettes, representing occluder surfaces as seen from the sources. The silhouettes are used to compute ray intersection during visibility computation.
*/

class SilhouetteContainer
{
public:
    SilhouetteContainer()
    {
    }

    ~SilhouetteContainer()
    {
        for (auto s : mSilhouettes)
            delete s;
    }

    const std::unordered_set<Silhouette*>& getSilhouettes()
    {
        return mSilhouettes;
    }

    void addSilhouette(Silhouette* aSilhouette)
    {
        if (mSilhouettes.find(aSilhouette) == mSilhouettes.end())
            mSilhouettes.insert(aSilhouette);
    }

    virtual bool intersect(VisibilityRay* aRay, double aDistance = 0)
    {
        GeometryRay myGeometryRay(*aRay);

        bool hasIntersection = false;
        for (auto s : mSilhouettes)
        {
            const auto& myMeshFaces = s->getMeshFaces();
            for (auto faceIndex : s->getSilhouetteFaces())
            {
                const SilhouetteMeshFace& face = myMeshFaces[faceIndex];
                bool hit = false;
                if (aDistance == 0.0)
                {
                    hit = MathGeometry::hitsTriangle<float>(myGeometryRay, face.getVertex(0), face.getVertex(1), face.getVertex(2));
                }
                else
                {
                    hit = MathGeometry::hitsCylinder<float>(myGeometryRay, aDistance, face.getVertex(0), face.getVertex(1), face.getVertex(2));
                }

                if (hit)
                {
                    aRay->addIntersection(s->getGeometryId(), faceIndex, 0.0);
                    hasIntersection = true;
                    if (aDistance == 0)
                    {
                        break;
                    }
                }
            }

        }

        return hasIntersection;
    }

    template<class P, class S>
    static bool isOccluded(PluckerPolytope<P>* polytope, PluckerPolyhedron<P>* polyhedron, const std::vector <Silhouette*>& aSilhouettes, const std::vector<P>& polytopeLines, S myTolerance)
    {
        for (auto iter = aSilhouettes.begin(); iter != aSilhouettes.end(); iter++)
        {
            Silhouette* s = (*iter);
            const auto& edgesProcessed = s->getEdgesProcessed();

            if (s->getAvailableEdgeCount() == 0)
            {
                bool allAreBlocked = true;

                for (size_t edgeIndex : edgesProcessed)
                {
                    SilhouetteEdge silhouetteEdge = s->getEdge(edgeIndex);
                    size_t hyperplaneIndex = silhouetteEdge.mHyperPlaneIndex;
                    V_ASSERT(hyperplaneIndex > 0);

                    const P& myHyperplane = polyhedron->get(hyperplaneIndex);

                    for (auto line : polytopeLines)
                    {
                        GeometryPositionType position = MathPredicates::getVertexPlaneRelativePosition(myHyperplane, line, myTolerance);
                        if (position != ON_NEGATIVE_SIDE)
                        {
                            allAreBlocked = false;
                            break;
                        }
                    }

                    if (!allAreBlocked)
                        break;
                }
                if (allAreBlocked)
                    return true;
            }
        }
        return false;
    }

    virtual void prepare() {};
private:
    std::unordered_set<Silhouette*> mSilhouettes;
};
}

