#pragma once
#include <random>
#include "RayTracer.h"
#include "ParticleMap.h"
#include "kdtree.h"

using namespace std;

#define RAYTRACE 0
#define PATHTRACE 1
#define PATHGUIDE 2

class Renderer {
public:
	Renderer (int numRays, int mode, RayTracer rayTracer) {
        m_numRays = numRays;
        m_mode = mode;
        m_rayTracer = rayTracer;
    }
	virtual ~Renderer() {}
    
    inline void shade (const Scene & scene,
                        const Ray & ray,
                        Vec3i & triangle,
                        size_t & meshIndex,
                        float & u,
                        float & v,
                        float & d,
                        Vec3f & hitNormal,
                        Vec3f & trianglePoint,
                        Vec3f & colorResponse) {
       float w = 1.f - u - v;
       const auto& mesh = scene.meshes()[meshIndex];
       const auto& P = mesh.vertexPositions();
       const auto& N = mesh.vertexNormals();
       const Camera& camera = scene.camera();
       //const Vec3i& triangle = mesh.indexedTriangles()[triangleIndex];
        hitNormal = normalize(w * N[triangle[0]] + u * N[triangle[1]] + v * N[triangle[2]]);
        trianglePoint = w * P[triangle[0]] + u * P[triangle[1]] + v * P[triangle[2]];
       Vec3f radiance, bsdf;
        
       Material material = mesh.material();
       Vec3f pointCameraDirection = camera.position() - trianglePoint;
        
        if (m_mode==PATHGUIDE){
            vector<Particle> result;
            int k = 200;
            Particle triangleParticle;
            triangleParticle.position() = trianglePoint;
            //m_photonKdTree.knearest(triangleParticle, k, result);
            // radius of the sphere containing the k photons
            float r = dist(result[k-1].position(), trianglePoint),
            area = M_PI * r * r;
            
            Vec3f averageDirection(0.f,0.f,0.f);
            for (Particle photon : result){
                averageDirection+=photon.weight() * photon.incomeDirection();
                radiance+=photon.weight()*Vec3f(1.f,1.f,1.f);
            }
            radiance /= area;
            bsdf = material.evaluateColorResponse(hitNormal, normalize(averageDirection), -ray.direction());
            colorResponse += radiance * bsdf;
        }
        else {
            for (LightSource lightSource : scene.lightsources()) {
               // Check that the point on the triangle is iluminated by light
               // Cast a ray from a point on the triangle to the light source
               Vec3f pointLightDirection = lightSource.randAreaPosition() - trianglePoint;
               Ray lightRay (trianglePoint,  pointLightDirection);
               if (m_rayTracer.rayTrace (lightRay, scene, meshIndex, triangle, u, v, d))
                   continue;
               bsdf = material.evaluateColorResponse(hitNormal, pointLightDirection, -ray.direction());
               radiance = lightSource.evaluateLight(trianglePoint);
               colorResponse += radiance * bsdf;
           }
        }
    }
    
	inline Vec3f calculateColorRay (const Scene & scene, Ray ray, bool & posIntersectionFound) {
        size_t meshIndex;
        Vec3i triangle;
        Vec3f colorResponse(0.f,0.f,0.f), hitNormal, trianglePoint;
       
        float u, v, d;
        bool intersectionFound = m_rayTracer.rayTrace (ray, scene, meshIndex, triangle, u, v, d);
        
        if (not (intersectionFound && d > 0.f)){
            posIntersectionFound = false;
            return colorResponse;
        }
        shade (scene, ray, triangle, meshIndex,  u,  v, d, hitNormal, trianglePoint, colorResponse);
        // check that all the values are less than one
        //for (int i = 0; i < 3; i++) colorResponse[i] = fmax(fmin(colorResponse[i], 1.f),0.f);
        return colorResponse;
	}
    

    inline Vec3f calculateColorPath (const Scene & scene, Ray ray, bool & posIntersectionFound, int depth, const int finalDepth, Vec3i sampleIds) {
        size_t meshIndex;
        Vec3i triangle;
        Vec3f colorResponse(0.f,0.f,0.f),  hitNormal, trianglePoint;
        if (depth >= finalDepth)
            return colorResponse;
        
        float u, v, d;
        bool intersectionFound = m_rayTracer.rayTrace (ray, scene, meshIndex, triangle, u, v, d);
        
        if (not (intersectionFound && d > 0.f)){
            if (depth==0) posIntersectionFound = false;
            return colorResponse;
        }
        shade (scene, ray, triangle, meshIndex,  u,  v, d, hitNormal, trianglePoint, colorResponse);
        
        // check that all the values are less than one
        
        //uniform_real_distribution<float> dis(-1.f,1.f);
        //Vec3f randomDirection(stratifiedSample1D(sampleIds[0], m_numRays, -1.f, 1.f),
                              //stratifiedSample1D(sampleIds[1], m_numRays, -1.f, 1.f),
                              //stratifiedSample1D(sampleIds[2], m_numRays, -1.f, 1.f));
        Vec3f randomDirection = m_rayTracer.hsphereUniformSample(hitNormal, M_PI / 2.f);
        Ray nextRay(trianglePoint, randomDirection);
        
        return colorResponse + calculateColorPath(scene, nextRay, posIntersectionFound, depth+1, finalDepth, sampleIds);
    }
    


	inline void render (const Scene& scene, Image& image) {
		size_t w = image.width();
		size_t h = image.height();
		const Camera& camera = scene.camera();
        if (m_mode==PATHGUIDE){
        // photon map test
         ParticleMap m_photonMap(scene, 50000, m_rayTracer);
         cout << m_photonMap.size() << endl;
        
         kdtree m_photonKdTree(m_photonMap.list.begin(), m_photonMap.list.end());
         int k = 200;
         
         Particle photon;
         photon.position() = { 0.2, 0, 0.426};
         vector<Particle> result;
         
         
         Particle n = m_photonKdTree.nearest(photon);
         
         std::cout << "Test kdtree:\n";
         std::cout << "nearest point: " << n.position() << '\n';
         std::cout << "distance: " << m_photonKdTree.distance() << '\n';
         m_photonKdTree.knearest(photon, k, result);
         for (int i = 0; i < k; i++)
              std::cout << "knearest point: " << i << ": " << result[i].position() << ", dist: " <<
              dist(photon.position(), result[i].position()) << endl;
         std::cout << "distance: " << m_photonKdTree.distance() << '\n';
         
         std::cout << "nodes visited: " << m_photonKdTree.visited() << '\n';

         
         m_photonMap.saveToPCD("pointcloud.pcd");
        }
        ///
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
                    Vec3f noise = m_rayTracer.jitterSample(i, m_numRays);
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
                        case PATHGUIDE:
                            colorResponse += normalizeColor(calculateColorRay (scene, ray, posIntersectionFound));
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
    ParticleMap m_photonMap,
    m_importonMap;
    RayTracer m_rayTracer;
    
    
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
