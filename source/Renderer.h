#pragma once
#include <random>
#include "RayTracer.h"
#include "PhotonMap.h"

#include "kdtree.h"

using namespace std;

#define RAYTRACE 0
#define PATHTRACE 1

class Renderer {
public:
    Renderer () {}
    
	Renderer (int numRays, int mode, RayTracer rayTracer) {
        m_numRays = numRays;
        m_mode = mode;
        m_rayTracer = rayTracer;
    }
    
    Renderer (int numRays, int mode, RayTracer rayTracer, int numPhotons, int k) {
        m_numRays = numRays;
        m_mode = mode;
        m_rayTracer = rayTracer;
        m_numPhotons = numPhotons;
        m_k = k;
        
    }
    
	virtual ~Renderer() {}
    
    inline Vec3f dotArr(const vector<Vec3f> & arr, Vec3i & triangle, float w, float u, float v){
        return w *arr[triangle[0]] + u * arr[triangle[1]] + v * arr[triangle[2]];
    }
    
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
        hitNormal = normalize(dotArr(N, triangle, w, u, v));
        trianglePoint = dotArr(P, triangle, w, u, v);
       Vec3f radiance, bsdf;
        
       Material material = mesh.material();
       Vec3f pointCameraDirection = camera.position() - trianglePoint;
        //cout <<"start"<< endl;
        
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
    
    inline void shade (const Scene & scene,
                     const Ray & ray,
                     Vec3i & triangle,
                     size_t & meshIndex,
                     float & u,
                     float & v,
                     float & d,
                     Vec3f & hitNormal,
                     Vec3f & trianglePoint,
                     Vec3f & colorResponse,
                     kdtree & particleTree) {
    float w = 1.f - u - v;
    const auto& mesh = scene.meshes()[meshIndex];
    const auto& P = mesh.vertexPositions();
    const auto& N = mesh.vertexNormals();
    const Camera& camera = scene.camera();
    //const Vec3i& triangle = mesh.indexedTriangles()[triangleIndex];
     hitNormal = normalize(dotArr(N, triangle, w, u, v));
     trianglePoint = dotArr(P, triangle, w, u, v);
    Vec3f radiance, bsdf;
     
    Material material = mesh.material();
    Vec3f pointCameraDirection = camera.position() - trianglePoint;
     //cout <<"start"<< endl;

        vector<Particle> result;
        Particle triangleParticle;
        triangleParticle.position() = trianglePoint;
        //cout <<"ooop"<< endl;
        particleTree.knearest(triangleParticle, m_k, result);
        //std::cout << "Test kdtree:\n";
        //for (int i = 0; i < m_k; i++)
             //std::cout << "knearest point: " << i << ": " << result[i].position() << ", dist: " <<
             //dist(triangleParticle.position(), result[i].position()) << endl;
        
        // radius of the sphere containing the k photons
        float r = dist(result[m_k-1].position(), trianglePoint),
        area = M_PI * r * r;
         
        Vec3f averageDirection(0.f,0.f,0.f);
        for (Particle photon : result){
            averageDirection+=photon.weight() * photon.incomeDirection();
            radiance+=photon.weight()*Vec3f(1.f,1.f,1.f);
        }
        //cout << radiance << endl;
        radiance /= area;
        radiance /= 7500;
        bsdf = material.evaluateColorResponse(hitNormal, normalize(averageDirection), -ray.direction());
        colorResponse += radiance * bsdf;
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
    
    inline Vec3f calculateColorRay (const Scene & scene, Ray ray, bool & posIntersectionFound, kdtree & particleTree) {
        size_t meshIndex;
        Vec3i triangle;
        Vec3f colorResponse(0.f,0.f,0.f), hitNormal, trianglePoint;
       
        float u, v, d;
        bool intersectionFound = m_rayTracer.rayTrace (ray, scene, meshIndex, triangle, u, v, d);
        
        if (not (intersectionFound && d > 0.f)){
            posIntersectionFound = false;
            return colorResponse;
        }
        shade (scene, ray, triangle, meshIndex,  u,  v, d, hitNormal, trianglePoint, colorResponse, particleTree);
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
    
    
    inline Vec3f calculateColorPath (const Scene & scene, Ray ray, bool & posIntersectionFound, int depth, const int finalDepth, Vec3i sampleIds, kdtree & particleTree) {
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
        shade (scene, ray, triangle, meshIndex,  u,  v, d, hitNormal, trianglePoint, colorResponse, particleTree);
        
        // check that all the values are less than one
        
        //uniform_real_distribution<float> dis(-1.f,1.f);
        //Vec3f randomDirection(stratifiedSample1D(sampleIds[0], m_numRays, -1.f, 1.f),
                              //stratifiedSample1D(sampleIds[1], m_numRays, -1.f, 1.f),
                              //stratifiedSample1D(sampleIds[2], m_numRays, -1.f, 1.f));
        Vec3f randomDirection = m_rayTracer.hsphereUniformSample(hitNormal, M_PI / 2.f);
        Ray nextRay(trianglePoint, randomDirection);
        
        return colorResponse + calculateColorPath(scene, nextRay, posIntersectionFound, depth+1, finalDepth, sampleIds, particleTree);
    }
    


	inline void render (const Scene& scene, Image& image) {
		size_t w = image.width();
		size_t h = image.height();
		const Camera& camera = scene.camera();
        // photon map test
            
         PhotonMap m_photonMap(scene, m_numPhotons, m_rayTracer);
        
         kdtree photonTree(m_photonMap.list().begin(), m_photonMap.list().end());
         //ParticleTree m_photonTree(m_photonMap);
         //cout << m_photonMap.size() << endl;
        
         //int k = 200;
         //vector<Particle> result;
            
         /*
          Particle photon;
         photon.position() = { 0.2, 0, 0.426};
         Particle n = m_photonMap.nearest(photon);
         
         std::cout << "Test kdtree:\n";
         std::cout << "nearest point: " << n.position() << '\n';
         std::cout << "distance: " << m_photonMap.tree().distance() << '\n';
         m_photonMap.knearest(photon, k, result);
         for (int i = 0; i < k; i++)
              std::cout << "knearest point: " << i << ": " << result[i].position() << ", dist: " <<
              dist(photon.position(), result[i].position()) << endl;
         std::cout << "distance: " << m_photonMap.tree().distance() << '\n';
         
         std::cout << "nodes visited: " << m_photonMap.tree().visited() << '\n';
        */
         
         //m_photonMap.saveToPCD("pointcloud.pcd");
        
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
                    Vec3f color(0.f, 0.f, 0.f);
                    switch(m_mode){
                        case RAYTRACE:
                            if (not photonTree.empty())
                                color = calculateColorRay (scene, ray, posIntersectionFound, photonTree);
                            else
                                color = calculateColorRay (scene, ray, posIntersectionFound);
                            break;
                        case PATHTRACE:
                            if (not photonTree.empty())
                                color = calculateColorPath (scene, ray, posIntersectionFound, 0, 3, Vec3i(vecInds[0][i], vecInds[1][i], vecInds[2][i]), photonTree);
                            else
                                color = calculateColorPath (scene, ray, posIntersectionFound, 0, 3, Vec3i(vecInds[0][i], vecInds[1][i], vecInds[2][i]));
                            break;
                        default:
                            break;
                    }
                    colorResponse += normalizeColor(color);
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
    int m_numRays, m_mode, m_numPhotons, m_k;
    PhotonMap m_photonMap,
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
