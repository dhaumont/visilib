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

#ifdef EMBREE

#include "embree3/rtcore.h"
#include "embree3/rtcore_ray.h"
#include "helper_triangle_mesh_container.h"

namespace visilib
{
    /** @brief Represents a set of silhouettes, using Embree library from Intel to perform ray-tracing against those silhouettes.

    This class is performance oriented, and is intended to provide the fastest possible occluder ray intersection, by using a modern ray tracing engine.
    */

    class SilhouetteContainerEmbree : public SilhouetteContainer
    {
    public:

        SilhouetteContainerEmbree();

        virtual ~SilhouetteContainerEmbree();

        /** @brief Intersect the ray with the geometry

        The function calls embree to perform the intersection
        */
        virtual bool intersect(VisibilityRay* aRay) override;


        /** @brief Prepare the scene before ray tracing

        The function commits the scene geometry to embree
        */
        virtual void prepare();
    private:

        RTCScene mScene;            /**< @brief The Embree scene*/
        static RTCDevice mDevice;   /**< @brief The Embree device*/
        std::unordered_map<unsigned int, std::pair<size_t, size_t>> ids;
    };

    inline SilhouetteContainerEmbree::SilhouetteContainerEmbree()
        :SilhouetteContainer()
    {
        if (mDevice == nullptr)
        {
            mDevice = rtcNewDevice("verbose=1");
           // rtcSetDeviceProperty(mDevice, RTC_DEVICE_PROPERTY_BACKFACE_CULLING_ENABLED, 1);
            //error_handler(rtcDeviceGetError(mDevice));

            /* set error handler */
            //	rtcDeviceSetErrorFunction(mDevice, error_handler);
        }
        mScene = rtcNewScene(mDevice);
    //    rtcSetSceneFlags(mScene, RTC_SCENE_FLAG_ROBUST);
        //rtcSetSceneFlags(mScene, RTC_CONFIG_BACKFACE_CULLING);
        rtcSetSceneFlags(mScene, RTC_SCENE_FLAG_ROBUST);
        rtcSetSceneFlags(mScene, RTC_SCENE_FLAG_CONTEXT_FILTER_FUNCTION);
    }

    inline SilhouetteContainerEmbree::~SilhouetteContainerEmbree()
    {
        rtcReleaseScene(mScene);
    }


#define HIT_LIST_LENGTH 16

    /** @brief Extended ray structure to store all the hits along the ray */
    struct GeometryRayMultiHit
    {
        RTCRayHit rayHit;

        // we remember up to 16 hits to ignore duplicate hits
        unsigned int firstHit, lastHit;
        unsigned int hit_geomIDs[HIT_LIST_LENGTH];
    };

    struct GeometryIntersectContextMultiHit
    {
        RTCIntersectContext context;
        GeometryRayMultiHit* ray;
    };

    inline RTCRayHit* RTCRayHit_(GeometryRayMultiHit& ray)
    {
        RTCRayHit* ray_ptr = (RTCRayHit*)& ray;
        return ray_ptr;
    }

    inline RTCRay* RTCRay_(GeometryRayMultiHit& ray)
    {
        RTCRay* ray_ptr = (RTCRay*)& ray;
        return ray_ptr;
    }


    inline void occlusionFilter(const RTCFilterFunctionNArguments* args)
    {
        /* avoid crashing when debug visualizations are used */
        if (args->context == nullptr) return;

        assert(args->N == 1);
        int* valid = args->valid;
        RTCHit * hit = (RTCHit*)args->hit;

        /* ignore inactive rays */
        if (valid[0] != -1) return;

        GeometryIntersectContextMultiHit* context = (GeometryIntersectContextMultiHit*)args->context;
        GeometryRayMultiHit* ray = context->ray;

        for (unsigned int i = ray->firstHit; i < ray->lastHit; i++)
        {
            unsigned slot = i % HIT_LIST_LENGTH;
            if (ray->hit_geomIDs[slot] == hit->geomID) {
                valid[0] = 0; return; // ignore duplicate intersections
            }
        }
        /* store hit in hit list */
        unsigned int slot = ray->lastHit % HIT_LIST_LENGTH;
        ray->hit_geomIDs[slot] = hit->geomID;
        ray->lastHit++;
        if (ray->lastHit - ray->firstHit >= HIT_LIST_LENGTH)
            ray->firstHit++;

        valid[0] = 0;
  }
    inline bool SilhouetteContainerEmbree::intersect(VisibilityRay* aRay)
    {
        GeometryRayMultiHit ray;
        ray.firstHit = 0;
        ray.lastHit = 0;
        ray.rayHit.ray.org_x = aRay->org[0];  ray.rayHit.ray.org_y = aRay->org[1]; ray.rayHit.ray.org_z = aRay->org[2];
        ray.rayHit.ray.dir_x = aRay->dir[0];  ray.rayHit.ray.dir_y = aRay->dir[1]; ray.rayHit.ray.dir_z = aRay->dir[2];
        ray.rayHit.ray.tnear = aRay->tnear;
        ray.rayHit.ray.tfar = aRay->tfar;
        ray.rayHit.ray.mask = -1;
        ray.rayHit.ray.time = 0;
        ray.rayHit.ray.id = 0;

        ray.rayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
        GeometryIntersectContextMultiHit context;
        rtcInitIntersectContext(&context.context);

        while (true)
        {
            context.ray = &ray;
            rtcIntersect1(mScene, &context.context, RTCRayHit_(ray));
            if (ray.rayHit.hit.geomID == RTC_INVALID_GEOMETRY_ID)
                break;
        }
        for (unsigned int i = ray.firstHit; i < ray.lastHit; i++)
        {
            auto h = ids[ray.hit_geomIDs[i]];
            aRay->addIntersection(h.first, h.second, 0.0);
        }
        return ray.lastHit > ray.firstHit;
     }

    inline void SilhouetteContainerEmbree::prepare()
    {
        for (auto s:getSilhouettes())
        {
            size_t id = s->getGeometryId();
            const auto& myMeshFaces = s->getMeshFaces();
            for (size_t faceIndex : s->getSilhouetteFaces())
            {
                RTCGeometry mesh = rtcNewGeometry(mDevice, RTC_GEOMETRY_TYPE_TRIANGLE);

                MathVector3f* vertices = (MathVector3f*)rtcSetNewGeometryBuffer(mesh, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, sizeof(MathVector3f), 3);
                for (int i = 0; i < 3; i++)
                {
                    vertices[i] = myMeshFaces[faceIndex].getVertex(i);
                }
#if 0
                MathVector3f center(vertices[0]);
                center += vertices[1];
                center += vertices[2];
                center /= 3.0;
                float s = 0.01;
                for (int i = 0; i < 3;i++)
                {
                    MathVector3f dx(vertices[i]); dx -= center; dx.normalize(); dx *= s; vertices[i] -= dx;
                }
#endif
                int* triangles = (int*)rtcSetNewGeometryBuffer(mesh, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, sizeof(int)*3, 1);
                triangles[0] = 0;
                triangles[1] = 1;
                triangles[2] = 2;
        //        triangles[3] = 0;
          //      triangles[4] = 2;
            //    triangles[5] = 1;

                rtcSetGeometryIntersectFilterFunction(mesh, occlusionFilter);

                rtcCommitGeometry(mesh);
                unsigned int geomID = rtcAttachGeometry(mScene, mesh);
                ids[geomID] = std::pair<size_t, size_t>(id, faceIndex);
                rtcReleaseGeometry(mesh);
            }
        }

        rtcCommitScene(mScene);
    }

}
#endif
