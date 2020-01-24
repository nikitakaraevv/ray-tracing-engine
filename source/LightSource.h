#pragma once
#include "Vec3.h"

class LightSource {
public:
	inline LightSource () {}
    inline LightSource (Vec3f position, Vec3f color, float intensity) {
        m_position = position;
        m_color = color;
        m_intensity = intensity;
    }
    
	virtual ~LightSource() {}

	inline Vec3f& position()  { return m_position; }

	inline Vec3f& color() { return m_color; }
	
    inline float intensity() { return m_intensity; }
    
    inline Vec3f evaluateLight (const Vec3f &point) {
        float d = dist(point, m_position);
        return m_color * (m_intensity / (ac + al * d + aq *d*d));
    }
private:
	Vec3f m_position,
          m_color;
    float m_intensity, ac = 0, al = 0, aq = 1;
    
};
