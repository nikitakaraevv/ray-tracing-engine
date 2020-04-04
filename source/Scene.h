#pragma once

#include <vector>

#include "Camera.h"
#include "LightSource.h"
#include "Mesh.h"

class Scene {
 public:
  inline Scene() {}
  virtual ~Scene() {}

  inline const Camera& camera() const { return m_camera; }

  inline Camera& camera() { return m_camera; }

  inline const std::vector<LightSource>& lightsources() const {
    return m_lightsources;
  }

  inline std::vector<LightSource>& lightsources() { return m_lightsources; }

  inline const std::vector<Mesh>& meshes() const { return m_meshes; }

  inline std::vector<Mesh>& meshes() { return m_meshes; }

 private:
  Camera m_camera;
  std::vector<LightSource> m_lightsources;
  std::vector<Mesh> m_meshes;
};
