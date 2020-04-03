#pragma once

#include <vector>
#include "Vec3.h"

using namespace std;

class Particle {
    public:
        Particle(){}
    
        inline Particle(Vec3f position, Vec3f direction, float weight) {
            m_position = position;
            m_direction = direction;
            m_weight = weight;
        }
    
        ~Particle(){}
    
        inline Vec3f& position() { return m_position; }
        
        inline const Vec3f& position() const { return m_position; }
    
        inline Vec3f& incomeDirection() { return m_direction; }
    
        inline const Vec3f& incomeDirection() const { return m_direction; }
    
        inline float& weight() { return m_weight; }
        
        inline const float& weight() const { return m_weight; }
    
    private:
        Vec3f m_position,
              m_direction;
        float m_weight;
};
