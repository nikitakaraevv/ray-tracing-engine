#pragma once

#include "Ray.h"
#include "Vec3.h"

/// Axis-aligned bounding box
class AABB {
 public:
  AABB(Vec3f minCorner = {0, 0, 0}, Vec3f maxCorner = {0, 0, 0})
      : m_minCorner(minCorner), m_maxCorner(maxCorner) {}

  inline const Vec3f& minCorner() const { return m_minCorner; }

  inline Vec3f& minCorner() { return m_minCorner; }

  inline const Vec3f& maxCorner() const { return m_maxCorner; }

  inline Vec3f& maxCorner() { return m_maxCorner; }

  bool isInsideTest(const Vec3f& v) const;

  bool intersectionTest(const Ray& r, Vec3f& entry, Vec3f& exit) const;

 private:
  Vec3f m_minCorner;
  Vec3f m_maxCorner;
};
