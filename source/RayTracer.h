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

using namespace std;

class RayTracer {
public:
	RayTracer () {}
	virtual ~RayTracer() {}
	   
	inline bool rayTrace (const Ray & ray, 
						  const Scene & scene, 
						  size_t & meshIndex, 
						  size_t & triangleIndex, 
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
				const Vec3i& triangle = T[tIndex];
				float ut, vt, dt;
				if (ray.triangleIntersect(P[triangle[0]], P[triangle[1]], P[triangle[2]], ut, vt, dt) == true) {
					if (dt > 0.f && dt < closest) {
						intersectionFound = true;
						closest = dt;
						meshIndex = mIndex;
						triangleIndex = tIndex;
						u = ut;
						v = vt;
						d = dt;
					}
				}
			}
		}
		return intersectionFound;
	}

	inline Vec3f shade (const Scene & scene, size_t meshIndex, size_t triangleIndex, float u, float v) {
		const auto& mesh = scene.meshes()[meshIndex];
		const auto& P = mesh.vertexPositions();
		const auto& N = mesh.vertexNormals();
        const Camera& camera = scene.camera();
		const Vec3i& triangle = mesh.indexedTriangles()[triangleIndex];
		Vec3f hitNormal = normalize((1.f - u - v) * N[triangle[0]] + u * N[triangle[1]] + v * N[triangle[2]]),
        trianglePoint = normalize((1.f - u - v) * P[triangle[0]] + u * P[triangle[1]] + v * P[triangle[2]]),
        radiance, bsdf;
        
        Material material = scene.material();
        LightSource lightSource = scene.lightsource();
        Vec3f shadow(0.f, 0.f, 0.f);
        
        // Check that the point on the triangle is iluminated by light
        // Cast a ray from a point on the triangle to the light source
        Vec3f pointLightDirection = lightSource.randAreaPosition() - trianglePoint,
              pointCameraDirection = camera.position() - trianglePoint;
        //cout << "position: " << lightSource.position() << endl;
        //cout << "randArea: " << lightSource.randAreaPosition() << endl;
        // if we are from the backside of the light source
        //if (dot(-pointLightDirection, lightSource.normal()) < 0.f)
        //    return shadow;
        Ray lightRay (trianglePoint,  pointLightDirection);
        float ut,vt,dt;
        size_t newTriangleIndex, newMeshIndex;
        // If the ray intersects other triangles, return shadow
        if (rayTrace (lightRay, scene, newMeshIndex, newTriangleIndex, ut, vt, dt))
            return shadow;
        
         
        
        bsdf = material.evaluateColorResponse(hitNormal, pointLightDirection, pointCameraDirection);
        radiance = lightSource.evaluateLight(trianglePoint);
        //cout <<"bsdf: " << bsdf << endl;
        cout <<"bsdf: " << bsdf << endl;
        cout <<"radiance: " <<  radiance << endl;
        cout <<"radiance*bsdf: " <<  radiance * bsdf << endl;
		
        return radiance * bsdf;//Vec3f(0.5f, 0.5f, 0.5f) + hitNormal / 2.f;
	}


	inline void render (const Scene& scene, Image& image) {
        
		size_t w = image.width();
		size_t h = image.height();
		const Camera& camera = scene.camera();
        
        uniform_real_distribution<float> dis(0.f,1.f);
        
		for (int y = 0; y < h; y++) {
			static int progress = 0;
			progress++;
			if (progress % 10 == 0)
				std::cout << ".";
#pragma omp parallel for
			for (int x = 0; x < w; x++) {
                // Trace N rays per pixel in ramdom directions to get rid of aliasing
                int N = 6, counter = 0;
                float shiftX, shiftY, eps = 1e-5;
                Vec3f colorResponse(0.f, 0.f, 0.f);
                for (int i = 0; i < N; i++){
                        shiftX = dis(gen);
                        shiftY = dis(gen);
                        Ray ray = camera.rayAt ((x + shiftX) / w, 1.f - (y + shiftY) / h);
                    
                        size_t meshIndex, triangleIndex;
                        float u, v, d;
                    bool intersectionFound = rayTrace (ray, scene, meshIndex, triangleIndex, u, v, d);
                    if (intersectionFound && d > 0.f){
                         colorResponse += shade (scene, meshIndex, triangleIndex, u, v);
                         counter++;
                    }
                }
                //cout << x << " " << y << endl;
                //cout << "colorResponse: " << colorResponse / float(N)<< endl;
                //cout << "image (x,y): "  << image (x,y) * ((N - counter)  / float(N)) << endl;
                image (x,y) = colorResponse  / float(N) +
                              image (x,y) * ((N - counter)  / float(N));
                
			}
		}
	}
    

    
        
};
