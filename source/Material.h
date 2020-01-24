#pragma once

#include "Vec3.h"
#include "math.h"

class Material {
public:
    inline Material () {
        m_kd = 1.3;
        m_albedo = Vec3f(0.5f, 0.5f, 0.5f);
    }
    inline Material (float kd, const Vec3f& albedo) {
        m_kd = kd;
        m_albedo = albedo;
    }
	virtual ~Material() {}

    inline float kd()  { return m_kd; }
    Vec3f evaluateColorResponse (const Vec3f& normal, const Vec3f& wi)  {
        float n_dot_wi = dot(normal, wi);
        return m_albedo * (m_kd / M_PI) * n_dot_wi;
    }
    

private:
    // kd is the diffuse coefficient
    float m_kd;
    Vec3f m_albedo;
};
