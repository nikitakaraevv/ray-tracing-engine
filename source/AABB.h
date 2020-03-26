#pragma once

#include "Vec3.h"
#include "Ray.h"

/// Axis-aligned bounding box
class AABB {
public:
    AABB(Vec3f minCorner={0,0,0}, Vec3f maxCorner={0,0,0}) {
            m_minCorner = minCorner;
            m_maxCorner = maxCorner;
        }
    
    inline const Vec3f& minCorner() const { return m_minCorner; }

    inline Vec3f& minCorner() { return m_minCorner; }
    
    inline const Vec3f& maxCorner() const { return m_maxCorner; }

    inline Vec3f& maxCorner() { return m_maxCorner; }
    
    bool isInsideTest(const Vec3f &v) const {
        for (int i = 0; i < 3; i++)
            if (v[i] < m_minCorner[i] - __FLT_EPSILON__ * 10 or
                v[i] > m_maxCorner[i] + __FLT_EPSILON__ * 10)
                return false;
        return true;
    }
    
    bool intersectionTest(const Ray &r, Vec3f &entry, Vec3f &exit) const {
        // Check if the ray origin is inside
        if (this->isInsideTest(r.origin())) entry = r.origin();
        else {
            // If the ray origin is not inside,
            // calculate intersections with the three candidate planes
            Vec3f t;
            for (int i = 0; i < 3; i++) {
                if (r.direction()[i] > 0)
                    t[i] = (m_minCorner[i] - r.origin()[i]) / r.direction()[i];
                else if (r.direction()[i] < 0)
                    t[i] = (m_maxCorner[i] - r.origin()[i]) / r.direction()[i];
                else
                    t[i] = -1;
            }
            float maxt = fmax(t[0], fmax(t[1], t[2]));
            if (maxt < 0) return false;
            
            // Check that intersection is inside AABB
            entry = r.origin() + maxt * r.direction();
            if (not this->isInsideTest(entry)) return false;
            
        }
        // So the ray intersects our AABB. We have also found an entry point.
        // Now we're looking for the exit point
        Vec3f t;
        for (int i = 0; i < 3; i++) {
            if (r.direction()[i] > 0)
                t[i] = (m_maxCorner[i] - r.origin()[i]) / r.direction()[i];
            else if (r.direction()[i] < 0)
                t[i] = (m_minCorner[i] - r.origin()[i]) / r.direction()[i];
            else
                t[i] = -1;
           
        }
        float mint = fmin(t[0], fmin(t[1], t[2]));
        exit = r.origin() + mint * r.direction();
        return true;
    }
    
private:
    Vec3f m_minCorner;
    Vec3f m_maxCorner;
};

