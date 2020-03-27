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

#define RAYTRACE 0
#define PATHTRACE 1
#define PATHGUIDE 2

class RayTracer {
public:
	RayTracer (int numRays, int mode) {
        m_numRays = numRays;
        m_mode = mode;
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
    

	inline Vec3f calculateColorRay (const Scene & scene, Ray ray, bool & posIntersectionFound) {
        size_t meshIndex;
        Vec3i triangle;
        Vec3f colorResponse(0.f,0.f,0.f);
        float w, u, v, d;
        bool intersectionFound = rayTrace (ray, scene, meshIndex, triangle, u, v, d);
        
        if (not (intersectionFound && d > 0.f)){
            posIntersectionFound = false;
            return colorResponse;
        }
        w = 1.f - u - v;
		const auto& mesh = scene.meshes()[meshIndex];
		const auto& P = mesh.vertexPositions();
		const auto& N = mesh.vertexNormals();
        const Camera& camera = scene.camera();
		//const Vec3i& triangle = mesh.indexedTriangles()[triangleIndex];
		Vec3f hitNormal = normalize(w * N[triangle[0]] + u * N[triangle[1]] + v * N[triangle[2]]),
        trianglePoint = normalize(w * P[triangle[0]] + u * P[triangle[1]] + v * P[triangle[2]]),
        radiance, bsdf;
        
        Material material = mesh.material();
        Vec3f pointCameraDirection = camera.position() - trianglePoint;
        for (LightSource lightSource : scene.lightsources()) {
            // Check that the point on the triangle is iluminated by light
            // Cast a ray from a point on the triangle to the light source
            Vec3f pointLightDirection = lightSource.randAreaPosition() - trianglePoint;
            Ray lightRay (trianglePoint,  pointLightDirection);
            if (rayTrace (lightRay, scene, meshIndex, triangle, u, v, d))
                continue;
            bsdf = material.evaluateColorResponse(hitNormal, pointLightDirection, pointCameraDirection);
            radiance = lightSource.evaluateLight(trianglePoint);
            colorResponse += radiance * bsdf;
        }
        // check that all the values are less than one
        //for (int i = 0; i < 3; i++) colorResponse[i] = fmax(fmin(colorResponse[i], 1.f),0.f);
        return colorResponse;
	}
    

    inline Vec3f calculateColorPath (const Scene & scene, Ray ray, bool & posIntersectionFound, int depth, const int finalDepth, Vec3i sampleIds) {
        size_t meshIndex;
        Vec3i triangle;
        Vec3f colorResponse(0.f,0.f,0.f);
        if (depth >= finalDepth)
            return colorResponse;
        
        float w, u, v, d;
        bool intersectionFound = rayTrace (ray, scene, meshIndex, triangle, u, v, d);
        
        if (not (intersectionFound && d > 0.f)){
            if (depth==0) posIntersectionFound = false;
            return colorResponse;
        }
        w = 1.f - u - v;
        const auto& mesh = scene.meshes()[meshIndex];
        const auto& P = mesh.vertexPositions();
        const auto& N = mesh.vertexNormals();
        const Camera& camera = scene.camera();
        //const Vec3i& triangle = mesh.indexedTriangles()[triangleIndex];
        Vec3f hitNormal = normalize(w * N[triangle[0]] + u * N[triangle[1]] + v * N[triangle[2]]),
        trianglePoint = normalize(w * P[triangle[0]] + u * P[triangle[1]] + v * P[triangle[2]]),
        radiance, bsdf;
        
        Material material = mesh.material();
        Vec3f pointCameraDirection = camera.position() - trianglePoint;
        for (LightSource lightSource : scene.lightsources()) {
            // Check that the point on the triangle is iluminated by light
            // Cast a ray from a point on the triangle to the light source
            Vec3f pointLightDirection = lightSource.randAreaPosition() - trianglePoint;
            Ray lightRay (trianglePoint,  pointLightDirection);
            if (rayTrace (lightRay, scene, meshIndex, triangle, u, v, d))
                continue;
            bsdf = material.evaluateColorResponse(hitNormal, pointLightDirection, pointCameraDirection);
            radiance = lightSource.evaluateLight(trianglePoint);
            colorResponse += radiance * bsdf;
        }
        // check that all the values are less than one
        
        uniform_real_distribution<float> dis(-1.f,1.f);
        //Vec3f randomDirection(stratifiedSample1D(sampleIds[0], m_numRays, -1.f, 1.f),
                              //stratifiedSample1D(sampleIds[1], m_numRays, -1.f, 1.f),
                              //stratifiedSample1D(sampleIds[2], m_numRays, -1.f, 1.f));
        Vec3f randomDirection(importanceSample(hitNormal));
        Ray nextRay(trianglePoint, randomDirection);
        
        return colorResponse + calculateColorPath(scene, nextRay, posIntersectionFound, depth+1, finalDepth, sampleIds);
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
                float shiftX, shiftY;
                Vec3f colorResponse(0.f, 0.f, 0.f);
                // create vector with indices for stratified sampling
                Vec3<vector<int>> vecInds;
                for (int i = 0; i < 3; i++){
                    vector<int> v(m_numRays); // vector with 100 ints.
                    std::iota (std::begin(v), std::end(v), 0);
                    vecInds[i] = v;
                }
                // shuffle them to get different indices for different coorfinates
                std::shuffle ( vecInds[0].begin(), vecInds[0].end(), gen );
                std::shuffle ( vecInds[1].begin(), vecInds[1].end(), gen );
                std::shuffle ( vecInds[2].begin(), vecInds[2].end(), gen );
                #pragma omp parallel for
                for (int i = 0; i < m_numRays; i++){
                    Vec3f noise = jitterSample(i, m_numRays);
                    shiftX = noise[0];
                    shiftY = noise[1];
                    Ray ray = camera.rayAt ((x + shiftX) / w, 1.f - (y + shiftY) / h);
                    bool posIntersectionFound = true;
                    switch(m_mode){
                        case RAYTRACE:
                            colorResponse += normalizeColor(calculateColorRay (scene, ray, posIntersectionFound));
                            break;
                        case PATHTRACE:
                            colorResponse += normalizeColor(calculateColorPath (scene, ray, posIntersectionFound, 0, 3, Vec3i(vecInds[0][i], vecInds[1][i], vecInds[2][i])));
                            break;
                        default:
                            break;
                    }
                    if (posIntersectionFound) counter++;
                }
                //cout << x << " " << y << endl;
                //cout << "colorResponse: " << colorResponse / float(N)<< endl;
                image (x,y) = colorResponse  / float(m_numRays) +
                              image (x,y) * ((m_numRays - counter)  / float(m_numRays));
                
			}
            
        }
    }
        
    
private:
    int m_numRays, m_mode;
    
    float stratifiedSample1D(int sampleIdx, int nSamples, float left, float right) {
        uniform_real_distribution<> dis(left, right);
        float invNSamples = (right-left) / nSamples,
        delta = 0.5f;
        return left + (sampleIdx + delta * dis(gen)) * invNSamples;
    }
    Vec3f importanceSample(Vec3f normal) {
        uniform_real_distribution<> dis(0.0, 1.0);
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

    
    inline Vec3f normalizeColor(Vec3f colorResponse) {
        for (int i = 0; i < 3; i++)
            colorResponse[i] = fmax(fmin(colorResponse[i], 1.f),0.f);
        return colorResponse;
    }
    
    void printProgressBar(float prop) {
        int progress = round(50.0f * prop);
        string progressBar = "";
        for (int i=0; i<progress; i++) {
            progressBar += "\u2588";
        }
        std::cout << "Raytracing... [" << progressBar << string(50 - progress, ' ') << "] " << progress * 2 << "%\r" << flush;
    }

};
