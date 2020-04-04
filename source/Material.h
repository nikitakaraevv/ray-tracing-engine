#pragma once

#include "Vec3.h"
#include "math.h"

class Material {
 public:
  inline Material() {
    m_kd = M_PI;
    m_albedo = Vec3f(0.9f, 0.4f, 0.4f);
    m_alpha = 0.5f;
    m_F0 = Vec3f(0.31f, 0.31f, 0.31f);  // glass;  gold: (1.0f, 0.86f, 0.57f)
  }
  inline Material(float kd, float alpha, const Vec3f& albedo, const Vec3f& F0) {
    // m_kd is the diffusion coefficient
    m_kd = kd;
    m_albedo = albedo;
    m_alpha = alpha;
    m_F0 = F0;
  }
  virtual ~Material() {}

  inline float kd() { return m_kd; }

  Vec3f evaluateColorResponse(const Vec3f& normal, const Vec3f& wi,
                              const Vec3f& wo) {
    Vec3f response =
        m_kd * diffuseResponse() +
        (1 - m_kd) *
            specularResponse(normalize(normal), normalize(wi), normalize(wo));

    for (int i = 0; i < 3; i++)
      if (response[i] < 0.f) response[i] = 0.f;

    return response;
  }

  Vec3f diffuseResponse() { return m_albedo / M_PI; }

  /// Microfacet specular reponse
  Vec3f specularResponse(const Vec3f& normal, const Vec3f& wi,
                         const Vec3f& wo) {
    Vec3f wh = wi + wo;
    wh = normalize(wh);
    // Distribution
    float D =
        m_alpha * m_alpha /
        (M_PI * pow(1 + (m_alpha * m_alpha - 1) * pow(dot(normal, wh), 2), 2));

    // Fresnel term
    Vec3f F =
        m_F0 + (Vec3f(1.f, 1.f, 1.f) - m_F0) * pow(1 - fmax(0, dot(wi, wh)), 5);
    // Geometric term
    float G = gSchlick(wi, normal) * gSchlick(wo, normal);

    float denom = (4. * dot(normal, wi) * dot(normal, wo));
    Vec3 microfacetTerm = D * F * G / denom;

    return microfacetTerm;
  }

 private:
  float m_kd, m_alpha;
  Vec3f m_albedo, m_F0;

  float gSchlick(const Vec3f& w, const Vec3f& normal) {
    float k = m_alpha * sqrt(2. / M_PI);
    return dot(normal, w) / (dot(normal, w) * (1 - k) + k);
  }
};
