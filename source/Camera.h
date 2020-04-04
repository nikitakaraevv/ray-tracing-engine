#pragma once

#include "Ray.h"
#include "Vec3.h"

/// Basic camera with "lookat" control
class Camera {
 public:
  Camera(const Vec3f& lookFrom = Vec3f(0.f, 0.f, 1.f),
         const Vec3f& lookAt = Vec3f(), const Vec3f& up = Vec3f(0.f, 1.f, 0.f),
         float verticalFoV = 60.f, float aspectRatio = 1.f)
      : m_verticalFoV(verticalFoV),
        m_aspectRatio(aspectRatio),
        m_position(lookFrom) {
    float angle = verticalFoV * M_PI / 180.f;
    float halfHeight = tan(angle / 2.f);
    float halfWidth = aspectRatio * halfHeight;
    Vec3f atFrom = normalize(lookFrom - lookAt);
    Vec3f u = normalize(cross(up, atFrom));
    Vec3f v = cross(atFrom, u);
    m_lowerLeftCorner = m_position - halfWidth * u - halfHeight * v - atFrom;
    m_horizontal = 2.f * halfWidth * u;
    m_vertical = 2.f * halfHeight * v;
  }

  /// Generate a ray for a given (u,v) coordinate on the image plane.
  Ray rayAt(float u, float v) const {
    return Ray(m_position, normalize(m_lowerLeftCorner + u * m_horizontal +
                                     v * m_vertical - m_position));
  }

  inline const Vec3f& position() const { return m_position; }

 private:
  float m_verticalFoV;
  float m_aspectRatio;
  Vec3f m_position;
  Vec3f m_lowerLeftCorner;
  Vec3f m_horizontal;
  Vec3f m_vertical;
};
