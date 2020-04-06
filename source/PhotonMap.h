#pragma once

#include <vector>

#include "RayTracer.h"
#include "Vec3.h"
#include "kdtree.h"

using namespace std;

class PhotonMap {
 public:
  PhotonMap() {}
  PhotonMap(const Scene& scene, int numOfPhotons, RayTracer rayTracer) {
    if (numOfPhotons > 0) {
      cout << "Constructing a photon map with " << numOfPhotons << " photons"
           << endl;
      m_rayTracer = rayTracer;
      float lightPdf = 1.f / scene.lightsources().size();
      int photonsPerLS = (int)(numOfPhotons * lightPdf);
      cout << "Emitting " << photonsPerLS << " photons per light source"
           << endl;

      // Save depths of emmited photons
      vector<int> depths(m_max_depth);
      for (LightSource lightSource : scene.lightsources()) {
        Vec3f lsNormal = lightSource.normal();
        for (int i = 0; i < photonsPerLS; i++) {
          //  Uniformly sample photons from an area light source
          Vec3f startPosition = lightSource.randAreaPosition(),
                startDirection =
                    m_rayTracer.hsphereUniformSample(lsNormal, M_PI / 2.f);

          //  Assign higher probabilities to photons with a direction close to
          //  normal
          float pdf = dot(normalize(startDirection), normalize(lsNormal)),
                weight = lightSource.radiance(startPosition) / (pdf * lightPdf);

          Ray ray(startPosition, startDirection);
          Particle photon({0., 0., 0.}, {0., 0., 0.}, weight);

          int depth = calculatePhotonPath(scene, ray, photon, 0, false);
          if (depth >= 0) depths[depth]++;
        }
      }
      for (int i = 0; i < depths.size(); i++)
        if (depths[i] > 0)
          cout << depths[i] << " photons with depth " << i << endl;
    }
  }

  int size() { return this->m_list.size(); }

  inline vector<Particle>& list() { return m_list; }

  inline const vector<Particle>& list() const { return m_list; }

  /// Save constructed photon map to a .pcd file (point cloud)
  inline void saveToPCD(const std::string& filename) {
    std::ofstream out(filename.c_str());
    if (!out) {
      std::cerr << "Cannot open file " << filename.c_str() << std::endl;
      std::exit(1);
    }
    out << "VERSION .7" << endl
        << "FIELDS x y z normal_x normal_y normal_z" << endl
        << "SIZE 4 4 4 4 4 4" << endl
        << "TYPE F F F F F F" << endl
        << "COUNT 1 1 1 1 1 1" << endl
        << "WIDTH " << m_list.size() << endl
        << "HEIGHT 1" << endl
        << "VIEWPOINT 0 0 0 1 0 0 0" << endl
        << "POINTS " << m_list.size() << endl
        << "DATA ascii" << endl;
    for (size_t i = 0; i < m_list.size(); i++) {
      Particle p(m_list[i]);
      out << p.position()[0] << " " << p.position()[1] << " " << p.position()[2]
          << " " << p.incomeDirection()[0] << " " << p.incomeDirection()[1]
          << " " << p.incomeDirection()[2] << " ";
      out << std::endl;
    }
    cout << "Particle map was saved to: " << filename << endl;
    out.close();
  }

 private:
  RayTracer m_rayTracer;
  vector<Particle> m_list;
  int m_max_depth = 20;

  /// Recursively emit photons
  inline int calculatePhotonPath(const Scene& scene, Ray ray, Particle photon,
                                 int depth, bool exit) {
    if (exit) {
      m_list.push_back(photon);
      return depth - 1;
    }
    if (depth >= m_max_depth) return -1;

    size_t meshIndex;
    Vec3i triangle;

    float w, u, v, d;
    bool intersectionFound =
        m_rayTracer.rayTrace(ray, scene, meshIndex, triangle, u, v, d);

    // save the last photon position and direction if an intersection is not
    // found
    if (not(intersectionFound && d > 0.f)) {
      if (depth != 0) m_list.push_back(photon);
      return -1;
    }
    w = 1.f - u - v;
    const auto& mesh = scene.meshes()[meshIndex];
    const auto& P = mesh.vertexPositions();
    const auto& N = mesh.vertexNormals();
    Material material = mesh.material();

    Vec3f hitNormal = normalize(w * N[triangle[0]] + u * N[triangle[1]] +
                                v * N[triangle[2]]),
          trianglePoint =
              w * P[triangle[0]] + u * P[triangle[1]] + v * P[triangle[2]];

    // update photon position and direction if an intersection is found
    photon.position() = trianglePoint;
    photon.incomeDirection() = -ray.direction();
    Vec3f randomDirection(
        m_rayTracer.hsphereUniformSample(hitNormal, M_PI / 2.f)),
        perfectReflection(ray.direction() -
                          2.f * (dot(ray.direction(), hitNormal)) * hitNormal);
    float bsdf = (material.evaluateColorResponse(hitNormal, ray.direction(),
                                                 randomDirection))
                     .length(),

          // pdf is the non-negative cos distance between perfect reflection and
          // a random direction
        pdf = (dot(normalize(randomDirection), normalize(perfectReflection)) +
               1.f) /
              2.f;
    photon.weight() *= bsdf / pdf;

    // russian roulette
    // probability depends on the weight of a photon
    float continueProb = fmin(photon.weight(), 1.f);
    uniform_real_distribution<float> dis(0.f, 1.f);

    if (dis(gen) > continueProb)
      exit = true;
    else
      photon.weight() /= continueProb;

    Ray nextRay(trianglePoint, randomDirection);

    return calculatePhotonPath(scene, nextRay, photon, depth + 1, exit);
  }
};
