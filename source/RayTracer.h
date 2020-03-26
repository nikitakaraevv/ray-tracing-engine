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

using namespace std;

void printProgressBar(float prop) {
    int progress = round(50.0f * prop);
    string progressBar = "";
    for (int i=0; i<progress; i++) {
        progressBar += "\u2588";
    }
    std::cout << "Raytracing... [" << progressBar << string(50 - progress, ' ') << "] " << progress * 2 << "%\r" << flush;
}

class RayTracer {
public:
	RayTracer (int numRays) {
        m_numRays = numRays;
    }
	virtual ~RayTracer() {}
    
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
            chrono::steady_clock::time_point begin = chrono::steady_clock::now();
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
            chrono::steady_clock::time_point end = chrono::steady_clock::now();
            int diff = chrono::duration_cast<chrono::nanoseconds>(end - begin).count();
            //if (diff>100000)
            //cout << "Raytr intersection took " << diff << "[ns]" << endl;
        }
        return intersectionFound;
    }
    
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
            chrono::steady_clock::time_point begin = chrono::steady_clock::now();
            vector<float> intersection = mesh.bvh().intersection(ray, triangle, mesh.vertexPositions(), mesh.indexedTriangles());
            chrono::steady_clock::time_point end = chrono::steady_clock::now();
            int diff = chrono::duration_cast<chrono::nanoseconds>(end - begin).count();
            //if (diff>100000)
            //cout << "Raytr intersection took " << diff << "[ns]" << endl;
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
    

	inline Vec3f shade (const Scene & scene, size_t meshIndex, Vec3i& triangle, float u, float v) {
		const auto& mesh = scene.meshes()[meshIndex];
		const auto& P = mesh.vertexPositions();
		const auto& N = mesh.vertexNormals();
        const Camera& camera = scene.camera();
		//const Vec3i& triangle = mesh.indexedTriangles()[triangleIndex];
		Vec3f hitNormal = normalize((1.f - u - v) * N[triangle[0]] + u * N[triangle[1]] + v * N[triangle[2]]),
        trianglePoint = normalize((1.f - u - v) * P[triangle[0]] + u * P[triangle[1]] + v * P[triangle[2]]),
        radiance, bsdf;
        Vec3f intersection_position(trianglePoint + 2 * __FLT_EPSILON__ * hitNormal);
        Material material = mesh.material();
        Vec3f colorResponse(0.f,0.f,0.f),
              pointCameraDirection = camera.position() - trianglePoint;
        for (LightSource lightSource : scene.lightsources()) {
            // Check that the point on the triangle is iluminated by light
            // Cast a ray from a point on the triangle to the light source
            Vec3f pointLightDirection = lightSource.randAreaPosition() - trianglePoint;
                  
            Ray lightRay (trianglePoint,  pointLightDirection);
            float ut,vt,dt;
            size_t  newMeshIndex;
            Vec3i newTriangle;
            
            if (rayTrace (lightRay, scene, newMeshIndex, newTriangle, ut, vt, dt))
                continue;
            bsdf = material.evaluateColorResponse(hitNormal, pointLightDirection, pointCameraDirection);
            radiance = lightSource.evaluateLight(trianglePoint);
            colorResponse += radiance * bsdf;

            //cout <<"bsdf: " << bsdf << endl;
            //cout <<"bsdf: " << bsdf << endl;
            //cout <<"radiance: " <<  radiance << endl;
            //cout <<"radiance*bsdf: " <<  radiance * bsdf << endl;
        }
        // check that all the values are less than one
        for (int i = 0; i < 3; i++) colorResponse[i] = fmax(fmin(colorResponse[i], 1.f),0.f);
        
        return colorResponse;//Vec3f(0.5f, 0.5f, 0.5f) + hitNormal / 2.f;
	}
    
    

	inline void render (const Scene& scene, Image& image) {
		size_t w = image.width();
		size_t h = image.height();
		const Camera& camera = scene.camera();
        
        uniform_real_distribution<float> dis(0.f,1.f);
		for (int y = 0; y < h; y++) {
            
			printProgressBar((float)(y + 1) / (float)w);
            #pragma omp parallel for
			for (int x = 0; x < w; x++) {
                // Trace m_numRays rays per pixel in random directions to get rid of aliasing
                int counter = 0;
                float shiftX, shiftY, eps = 1e-5;
                Vec3f colorResponse(0.f, 0.f, 0.f);
                //
                #pragma omp parallel for
                for (int i = 0; i < m_numRays; i++){
                    shiftX = dis(gen);
                    shiftY = dis(gen);
                    Ray ray = camera.rayAt ((x + shiftX) / w, 1.f - (y + shiftY) / h);
                
                    size_t meshIndex;
                    Vec3i triangle;
                    float u, v, d;
                    
                    bool intersectionFound = rayTrace (ray, scene, meshIndex, triangle, u, v, d);
                    
                    if (intersectionFound && d > 0.f){
                        colorResponse += shade (scene, meshIndex, triangle, u, v);
                        // colorResponse += pathShade (scene, meshIndex, triangle, u, v);
                        counter++;
                        
                    }
                }
                //cout << x << " " << y << endl;
                //cout << "colorResponse: " << colorResponse / float(N)<< endl;
                //cout << "image (x,y): "  << image (x,y) * ((N - counter)  / float(N)) << endl;
                image (x,y) = colorResponse  / float(m_numRays) +
                              image (x,y) * ((m_numRays - counter)  / float(m_numRays));
                
			}
            
        }
    }
        
    
private:
    int m_numRays;
};
