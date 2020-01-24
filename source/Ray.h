#pragma once

#include <limits>
#include "Vec3.h"

/// A 3D ray with basic intersection primitives for ray tracing.
class Ray {
public:
    inline Ray (const Vec3f & origin, const Vec3f & direction) : m_origin (origin), m_direction (direction) {}
    
    inline ~Ray () {}

    inline const Vec3f & origin () const { return m_origin; }
    
    inline const Vec3f & direction () const { return m_direction; }
    
	/// Return true if a the ray intersection the triangle [p0,p1,p2]. In this case, (u,v) are filled with barycentric coordinates (with the third one w = 1.0 - u - v), t with the distance from the origin to the intersection 
    bool triangleIntersect (const Vec3f &p0,
                            const Vec3f &p1,
                            const Vec3f &p2,
                            float & u,
                            float & v,
                            float & t) const;

private:
    Vec3f m_origin;
    Vec3f m_direction;
};

