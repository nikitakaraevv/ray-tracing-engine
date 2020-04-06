#pragma once
#include "Renderer.h"

#include <random>
#include <string>

#include "PhotonMap.h"
#include "RayTracer.h"
#include "kdtree.h"
using namespace std;

#define RAYTRACE 0
#define PATHTRACE 1

Renderer::Renderer(Scene& scene, int numRays, int mode, RayTracer rayTracer) {
  m_scene = scene;
  m_numRays = numRays;
  m_mode = mode;
  m_rayTracer = rayTracer;
  m_numPhotons = 0;
}

Renderer::Renderer(Scene& scene, int numRays, int mode, RayTracer rayTracer,
                   int numPhotons, int k) {
  m_scene = scene;
  m_numRays = numRays;
  m_mode = mode;
  m_rayTracer = rayTracer;
  m_numPhotons = numPhotons;
  m_k = k;
}

void Renderer::shade(const Ray& ray, Vec3i& triangle, size_t& meshIndex,
                     float& u, float& v, float& d, Vec3f& hitNormal,
                     Vec3f& trianglePoint, Vec3f& colorResponse) {
  float w = 1.f - u - v;
  const auto& mesh = m_scene.meshes()[meshIndex];
  const auto& P = mesh.vertexPositions();
  const auto& N = mesh.vertexNormals();
  const Camera& camera = m_scene.camera();

  hitNormal = normalize(dotArr(N, triangle, w, u, v));
  trianglePoint = dotArr(P, triangle, w, u, v);
  Vec3f radiance, bsdf;

  Material material = mesh.material();
  Vec3f pointCameraDirection = camera.position() - trianglePoint;

  for (LightSource lightSource : m_scene.lightsources()) {
    // Check that the point on the triangle is iluminated by light
    // Cast a ray from a point on the triangle to the light source
    Vec3f pointLightDirection = lightSource.randAreaPosition() - trianglePoint;
    Ray lightRay(trianglePoint, pointLightDirection);
    if (m_rayTracer.rayTrace(lightRay, m_scene, meshIndex, triangle, u, v, d))
      continue;
    bsdf = material.evaluateColorResponse(hitNormal, pointLightDirection,
                                          -ray.direction());
    radiance = lightSource.evaluateLight(trianglePoint);
    colorResponse += radiance * bsdf;
  }
}

void Renderer::shade(const Ray& ray, Vec3i& triangle, size_t& meshIndex,
                     float& u, float& v, float& d, Vec3f& hitNormal,
                     Vec3f& trianglePoint, Vec3f& colorResponse,
                     kdtree& particleTree) {
  float w = 1.f - u - v;
  const auto& mesh = m_scene.meshes()[meshIndex];
  const auto& P = mesh.vertexPositions();
  const auto& N = mesh.vertexNormals();
  const Camera& camera = m_scene.camera();

  hitNormal = normalize(dotArr(N, triangle, w, u, v));
  trianglePoint = dotArr(P, triangle, w, u, v);
  Vec3f radiance, bsdf;

  Material material = mesh.material();
  Vec3f pointCameraDirection = camera.position() - trianglePoint;

  vector<Particle> result;
  Particle triangleParticle;
  triangleParticle.position() = trianglePoint;

  // computing radiance using photon map
  particleTree.knearest(triangleParticle, m_k, result);

  // radius of the sphere containing the k photons
  float r = dist(result[m_k - 1].position(), trianglePoint),
        area = M_PI * r * r;

  Vec3f averageDirection(0.f, 0.f, 0.f);
  for (Particle photon : result) {
    averageDirection += photon.incomeDirection();  // photon.weight() *
                                                   // photon.incomeDirection();
    radiance +=
        Vec3f(1.f, 1.f, 1.f);  // photon.weight() * Vec3f(1.f, 1.f, 1.f);
  }
  radiance /= area;
  radiance /= m_numPhotons;
  radiance *= m_factor;
  bsdf = material.evaluateColorResponse(hitNormal, normalize(averageDirection),
                                        -ray.direction());
  colorResponse += radiance * bsdf;
}

Vec3f Renderer::calculateColorRay(Ray ray, bool& posIntersectionFound) {
  size_t meshIndex;
  Vec3i triangle;
  Vec3f colorResponse(0.f, 0.f, 0.f), hitNormal, trianglePoint;

  float u, v, d;
  bool intersectionFound =
      m_rayTracer.rayTrace(ray, m_scene, meshIndex, triangle, u, v, d);

  if (not(intersectionFound && d > 0.f)) {
    posIntersectionFound = false;
    return colorResponse;
  }
  shade(ray, triangle, meshIndex, u, v, d, hitNormal, trianglePoint,
        colorResponse);
  return colorResponse;
}

Vec3f Renderer::calculateColorRay(Ray ray, bool& posIntersectionFound,
                                  kdtree& particleTree) {
  size_t meshIndex;
  Vec3i triangle;
  Vec3f colorResponse(0.f, 0.f, 0.f), hitNormal, trianglePoint;

  float u, v, d;
  bool intersectionFound =
      m_rayTracer.rayTrace(ray, m_scene, meshIndex, triangle, u, v, d);

  if (not(intersectionFound && d > 0.f)) {
    posIntersectionFound = false;
    return colorResponse;
  }
  shade(ray, triangle, meshIndex, u, v, d, hitNormal, trianglePoint,
        colorResponse, particleTree);
  return colorResponse;
}

Vec3f Renderer::calculateColorPath(Ray ray, bool& posIntersectionFound,
                                   int depth, const int finalDepth) {
  size_t meshIndex;
  Vec3i triangle;
  Vec3f colorResponse(0.f, 0.f, 0.f), hitNormal, trianglePoint;
  if (depth >= finalDepth) return colorResponse / (float)depth;

  float u, v, d;
  bool intersectionFound =
      m_rayTracer.rayTrace(ray, m_scene, meshIndex, triangle, u, v, d);

  if (not(intersectionFound && d > 0.f)) {
    if (depth == 0)
      posIntersectionFound = false;
    else
      colorResponse /= (float)depth;
    return colorResponse;
  }
  shade(ray, triangle, meshIndex, u, v, d, hitNormal, trianglePoint,
        colorResponse);

  Vec3f randomDirection =
      m_rayTracer.hsphereUniformSample(hitNormal, M_PI / 2.f);
  Ray nextRay(trianglePoint, randomDirection);

  return colorResponse + calculateColorPath(nextRay, posIntersectionFound,
                                            depth + 1, finalDepth);
}

Vec3f Renderer::calculateColorPath(Ray ray, bool& posIntersectionFound,
                                   int depth, const int finalDepth,
                                   kdtree& particleTree) {
  size_t meshIndex;
  Vec3i triangle;
  Vec3f colorResponse(0.f, 0.f, 0.f), hitNormal, trianglePoint;
  if (depth >= finalDepth) return colorResponse / (float)depth;

  float u, v, d;
  bool intersectionFound =
      m_rayTracer.rayTrace(ray, m_scene, meshIndex, triangle, u, v, d);

  if (not(intersectionFound && d > 0.f)) {
    if (depth == 0)
      posIntersectionFound = false;
    else
      colorResponse /= (float)depth;
    return colorResponse;
  }
  shade(ray, triangle, meshIndex, u, v, d, hitNormal, trianglePoint,
        colorResponse, particleTree);

  Vec3f randomDirection =
      m_rayTracer.hsphereUniformSample(hitNormal, M_PI / 2.f);
  Ray nextRay(trianglePoint, randomDirection);

  return colorResponse + calculateColorPath(nextRay, posIntersectionFound,
                                            depth + 1, finalDepth,
                                            particleTree);
}

void Renderer::render(Image& image) {
  size_t w = image.width();
  size_t h = image.height();
  const Camera& camera = m_scene.camera();
  // photon map test
  Image updateImage(w, h), saveImage(w, h);
  PhotonMap m_photonMap(m_scene, m_numPhotons, m_rayTracer);
  m_photonMap.saveToPCD("pointcloud.pcd");
  if (m_numPhotons > 0)
    cout << "Constructing a kd-tree for the photon map." << endl;
  kdtree photonTree(m_photonMap.list().begin(), m_photonMap.list().end());

  vector<vector<int> > counter(w, vector<int>(h));

  uniform_real_distribution<float> dis(0.f, 1.f);
  // loop through number of rays
  for (int i = 0; i < m_numRays; i++) {
#pragma omp parallel for
    // loop through vertical image pixels
    for (int y = 0; y < h; y++) {
      printProgressBar((i * h + (float)(y + 1)) / ((float)h * m_numRays));

#pragma omp parallel for
      // loop through horizontal image pixels
      for (int x = 0; x < w; x++) {
        Vec3f colorResponse(0.f, 0.f, 0.f);
        Vec3f noise = m_rayTracer.jitterSample(i, m_numRays);

        float shiftX = noise[0];
        float shiftY = noise[1];
        Ray ray = camera.rayAt((x + shiftX) / (float)w,
                               1.f - (y + shiftY) / (float)h);
        bool posIntersectionFound = true;
        Vec3f color(0.f, 0.f, 0.f);
        switch (m_mode) {
          case RAYTRACE:
            if (not photonTree.empty())
              color = calculateColorRay(ray, posIntersectionFound, photonTree);
            else
              color = calculateColorRay(ray, posIntersectionFound);
            break;
          case PATHTRACE:
            if (not photonTree.empty())
              color = calculateColorPath(ray, posIntersectionFound, 0, 3,
                                         photonTree);
            else
              color = calculateColorPath(ray, posIntersectionFound, 0, 3);
            break;
          default:
            break;
        }
        colorResponse = normalizeColor(color);
        if (posIntersectionFound) {
          counter[x][y]++;
        }
        updateImage(x, y) += colorResponse;
      }
    }
    // Save an updated image after each iteration
    for (int y = 0; y < h; y++)
      for (int x = 0; x < w; x++)
        saveImage(x, y) = (updateImage(x, y) / float(i + 1)) +
                          image(x, y) * (i + 1 - counter[x][y]) / float(i + 1);
    // std::string saveName = "photonmapPathtrace/update_";
    // saveName += to_string(i) + ".ppm";
    string saveName = "update.ppm";
    saveImage.savePPM(saveName);
  }
  image = saveImage;
}

inline Vec3f Renderer::dotArr(const vector<Vec3f>& arr, Vec3i& triangle,
                              float w, float u, float v) {
  return w * arr[triangle[0]] + u * arr[triangle[1]] + v * arr[triangle[2]];
}

inline Vec3f Renderer::normalizeColor(Vec3f colorResponse) {
  for (int i = 0; i < 3; i++)
    colorResponse[i] = fmax(fmin(colorResponse[i], 1.f), 0.f);
  return colorResponse;
}

void Renderer::printProgressBar(float prop) {
  int progress = round(50.0f * prop);
  string progressBar = "";
  for (int i = 0; i < progress; i++) {
    progressBar += "\u2588";
  }
  std::cout << "Raytracing... [" << progressBar << string(50 - progress, ' ')
            << "] " << progress * 2 << "%\r" << flush;
}
