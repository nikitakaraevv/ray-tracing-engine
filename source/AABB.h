#pragma once

#include "Vec3.h"
#include "Ray.h"

/// Axis-aligned bounding box
class AABB {
public:
    AABB(Vec3f minCorner,
         Vec3f maxCorner)  {
        m_minCorner = minCorner;
        m_maxCorner = maxCorner;
    }
    bool isInsideTest(const Vec3f &v){
        for (int i = 0; i < 3; i++)
            if (v[i] < m_minCorner[i] or
                v[i] > m_maxCorner[i])
                return false;
        return true;
    }
    
    bool intersectionTest(const Ray &r, Vec3f &entry, Vec3f &exit) {
        // Check if the ray origin is inside
        if (isInsideTest(r.origin())) entry = r.origin();
        else {
            // If the ray origin is not inside,
            // calculate intersections with the three candidate planes
            Vec3f t;
            for (int i = 0; i < 3; i++) {
                if (r.direction()[i] > 0)
                    t[i] = (m_minCorner[i] - r.origin()[i]) / r.direction()[i];
                else
                    t[i] = (m_maxCorner[i] - r.origin()[i]) / r.direction()[i];
            }
            float maxt = max(t[0], max(t[1], t[2]));
            if (maxt < 0) return false;
            
            // Check that intersection is inside AABB
            entry = r.origin() + maxt * r.direction();
            if (not isInsideTest(entry)) return false;
            
        }
        // So the ray intersects our AABB. We have also found an entry point.
        // Now we're looking for the exit point
        Vec3f t;
        for (int i = 0; i < 3; i++) {
            if (r.direction()[i] > 0)
                t[i] = (m_maxCorner[i] - r.origin()[i]) / r.direction()[i];
            else
                t[i] = (m_minCorner[i] - r.origin()[i]) / r.direction()[i];
           
        }
        float mint = min(t[0], min(t[1], t[2]));
        exit = r.origin() + mint * r.direction();
        return true;
    }
    
private:
    Vec3f m_minCorner;
    Vec3f m_maxCorner;
};

