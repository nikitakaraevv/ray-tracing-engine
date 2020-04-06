#pragma once
#include <random>
#include <string>

#include "PhotonMap.h"
#include "RayTracer.h"
#include "kdtree.h"
using namespace std;

#define RAYTRACE 0
#define PATHTRACE 1

///  The most important class. All the rendering algorithms are implemented
///  here.
class Renderer {
 public:
  Renderer() {}
  ///  Constructor for naive ray-tracing and path-tracing algorithms
  Renderer(Scene& scene, int numRays, int mode, RayTracer rayTracer);

  ///  Constructor for photon mapping-based ray-tracing and path-tracing
  ///  algorithms
  Renderer(Scene& scene, int numRays, int mode, RayTracer rayTracer,
           int numPhotons, int k);

  virtual ~Renderer() {}

  /**
   *  Iteratively traces rays through image pixels.
   *  @param m_numRays - number of traced rays through pixel
   *  @param m_numPhotons - number of photons in the photon map.
   *  Photon map is built if  number of photons is greater than zero.
   *  @param m_k - number of neighbours in the photon map
   *  Updates rendered image after each ray traced through all the pixels
   */
  void render(Image& image);

  void savePhotonMap();

 private:
  int m_numRays, m_mode, m_numPhotons, m_k;
  PhotonMap m_photonMap, m_importonMap;
  RayTracer m_rayTracer;
  Scene m_scene;
  float m_factor = 100.f;

  ///  Shading for naive ray-tracing and path-tracing algorithms
  ///  Computes shades for each light source, computes radiance and BSDF.
  ///  Saves the color response
  void shade(const Ray& ray, Vec3i& triangle, size_t& meshIndex, float& u,
             float& v, float& d, Vec3f& hitNormal, Vec3f& trianglePoint,
             Vec3f& colorResponse);

  ///  Shading for photon mapping-based ray-tracing and path-tracing algorithms
  void shade(const Ray& ray, Vec3i& triangle, size_t& meshIndex, float& u,
             float& v, float& d, Vec3f& hitNormal, Vec3f& trianglePoint,
             Vec3f& colorResponse, kdtree& particleTree);

  ///  Calls shade to compute color response for a ray if an intersection is
  ///  found. Function is used for naive ray-tracing and path-tracing algorithms
  Vec3f calculateColorRay(Ray ray, bool& posIntersectionFound);

  ///  Function is used for photon mapping-based ray-tracing and path-tracing
  ///  algorithms
  Vec3f calculateColorRay(Ray ray, bool& posIntersectionFound,
                          kdtree& particleTree);

  ///  Recursively calls shade to compute color response for a path if an
  ///  intersection is found. Function is used for naive ray-tracing and
  ///  path-tracing algorithms
  Vec3f calculateColorPath(Ray ray, bool& posIntersectionFound, int depth,
                           const int finalDepth);

  ///  Function is used for photon mapping-based ray-tracing and path-tracing
  ///  algorithms
  Vec3f calculateColorPath(Ray ray, bool& posIntersectionFound, int depth,
                           const int finalDepth, kdtree& particleTree);

  void printProgressBar(float prop);

  Vec3f normalizeColor(Vec3f colorResponse);

  Vec3f dotArr(const vector<Vec3f>& arr, Vec3i& triangle, float w, float u,
               float v);
};
