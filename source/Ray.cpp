#include "Ray.h"

#include <algorithm>

using namespace std;

static float EPSILON = 0.000001f;

bool Ray::triangleIntersect (const Vec3f &p0,
                             const Vec3f &p1,
                             const Vec3f &p2,
                             float & u,
                             float & v,
                             float & t) const {
    Vec3f edge1 = p1 - p0, edge2 = p2 - p0;
    Vec3f pvec = cross(m_direction, edge2);
    float det = dot(edge1, pvec);
    if (fabs (det) < EPSILON)
        return false;
    float inv_det = 1.0f / det;
    Vec3f tvec = m_origin - p0;
    u = dot(tvec, pvec) * inv_det;
    Vec3f qvec = cross(tvec, edge1);
    v = dot(m_direction, qvec) * inv_det;
    t = dot(edge2, qvec) * inv_det;
    if (u < 0.f || u > 1.f)
        return false;
    if (v >= 0.f && u + v <= 1.f)
        return true;
    return false;
}