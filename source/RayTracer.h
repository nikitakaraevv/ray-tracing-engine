#pragma once
#include <random>
#include <cmath>
#include <algorithm>
#include <limits>
#include <random>
#include "Vec3.h"
#include "Image.h"
#include "Camera.h"
#include "Scene.h"
#include <chrono>
#include "Particle.h"

using namespace std;

///  This class allows to trace rays either in a brute force manner or using BVH
class RayTracer {
public:
	RayTracer () {
    }
	virtual ~RayTracer() {}
    
    /**
     * Traces a ray. If an intersection is found, returns true.
     * Saves index of the intersected mesh, intersected triangle,
     * and barycentric coordinates on it.
    */
    inline bool rayTrace (const Ray & ray,
                          const Scene & scene,
                          size_t & meshIndex,
                          Vec3i & triangle,
                          float & u,
                          float & v,
                          float & d) {
        const auto& meshes = scene.meshes();
        float closest = std::numeric_limits<float>::max();
        bool intersectionFound = false;
        for (size_t mIndex = 0; mIndex < meshes.size(); mIndex++) {
            const auto& P = meshes[mIndex].vertexPositions();
            const auto& T = meshes[mIndex].indexedTriangles();
            for (size_t tIndex = 0; tIndex < T.size(); tIndex++) {
                const Vec3i& temp_triangle = T[tIndex];
                float ut, vt, dt;
                if (ray.triangleIntersect(P[temp_triangle[0]], P[temp_triangle[1]], P[temp_triangle[2]], ut, vt, dt) == true) {
                    if (dt > 0.f && dt < closest) {
                        intersectionFound = true;
                        closest = dt;
                        triangle = temp_triangle;
                        meshIndex = mIndex;
                        u = ut;
                        v = vt;
                        d = dt;
                    }
                }
            }
        }
        return intersectionFound;
    }
    
    /// The same function using BVH
	inline bool rayTraceBVH (Ray & ray,
                             const Scene & scene,
                             size_t & nearest_index,
                             Vec3i & nearest_triangle,
                             float & u,
                             float & v,
                             float & d) {
        
        vector<Mesh> meshes = scene.meshes();
        int mesh_index = 0;
        vector<float> nearest_intersection = {};
        for (Mesh const &mesh : meshes) {
            Vec3i triangle;
            vector<float> intersection = mesh.bvh().intersection(ray, triangle, mesh.vertexPositions(), mesh.indexedTriangles());

            if (intersection.size() > 0) {
                if (nearest_intersection.size() == 0 or nearest_intersection[3] > intersection[3]) {
                        nearest_intersection = intersection;
                        nearest_triangle = triangle;
                        nearest_index = mesh_index;
                }
            }
            mesh_index++;
        }
        
        if (nearest_intersection.size() > 0 and nearest_intersection[3]>0.f) {
            u =  nearest_intersection[1];
            v =  nearest_intersection[2];
            d =  nearest_intersection[3];
            return true;
        }
        return false;
    }
    
    /// Different sampling functions
    float stratifiedSample1D(int sampleIdx, int nSamples, float left, float right) {
        uniform_real_distribution<> dis(left, right);
        float invNSamples = (right-left) / nSamples,
        delta = 0.5f;
        return left + (sampleIdx + delta * dis(gen)) * invNSamples;
    }
    
    Vec3f hsphereUniformSample(Vec3f normal, float maxRayAngle) {
        uniform_real_distribution<> dis(0.0, 2 * maxRayAngle / M_PI);
        normal.normalize();
        Vec3f v1, v2;
        normal.getTwoOrthogonals(v1, v2);
        v1.normalize();
        v2.normalize();
        float theta = asin(dis(gen));
        float phi = 2 * M_PI * dis(gen);
        Vec3f direction = v1 * cos(phi) + v2 * sin(phi);
        direction.normalize();
        return normalize(normal * cos(theta) + direction * sin(theta));
    }
    
    
    Vec3f jitterSample(int sampleIdx, int nSamples) {
        uniform_real_distribution<> dis(0.0, 1.0);
        int d = int(sqrt(float(nSamples)));
        int j2 = sampleIdx / d;
        int i2 = sampleIdx % d;
        float x = (float(i2) + dis(gen)) / float(d);
        float y = (float(j2) + dis(gen)) / float(d);
        return {x, y, 0};
    }

};
