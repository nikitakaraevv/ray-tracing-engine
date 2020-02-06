#pragma once
#include "Vec3.h"
#include <random>
#include <algorithm>
std::default_random_engine gen;

using namespace std;
class LightSource {
public:
	inline LightSource () {}
    inline LightSource (Vec3f position, Vec3f color, float intensity) {
        m_position = position;
        m_color = color;
        m_intensity = intensity;
    }
    inline LightSource (Vec3f position,
                        Vec3f color,
                        Vec3f direction,
                        float intensity,
                        float sideLength) : uniform(-sideLength, sideLength)
    {
        m_position = position;
        m_color = color;
        m_intensity = intensity;
        m_sideLength = sideLength;
        m_direction = direction;
        // calculate the normal and basis vectors for our square area
        m_normal = normalize(m_direction - m_position);
        m_vertical = normalize(cross(m_normal, m_normal + Vec3f(1.f,0.f,0.f)));
        m_horizontal = normalize(cross(m_normal, m_vertical));
        //cout << "m_vertical: " << m_vertical <<endl;
        //cout << "m_horizontal: " << m_horizontal <<endl;
    }
    
	virtual ~LightSource() {}

	inline Vec3f& position()  { return m_position; }
    
    inline Vec3f& normal()  { return m_normal; }
    
    inline Vec3f randAreaPosition() {
        return m_position + (uniform(gen) * m_vertical) + (uniform(gen) * m_horizontal);
    }

	inline Vec3f& color() { return m_color; }
	
    inline float intensity() { return m_intensity; }
    
    inline Vec3f evaluateLight (const Vec3f &point) {
        float d = dist(point, m_position);
        cout << "dist: " << d << endl;
        return 0.3f * m_color  + Vec3f(1.f, 1.f, 1.f) * fmin(0.7f, m_intensity / (ac + al * d + aq * d*d));
    }
private:
	Vec3f m_position,
          m_color,
          m_direction,
          m_normal,
          m_vertical,
          m_horizontal;
    
    uniform_real_distribution<float> uniform;
    float m_intensity, m_sideLength = 0., ac = 0, al = 0, aq = 1;
    
    
};
