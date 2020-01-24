#pragma once

#include <vector>

#include "Camera.h"
#include "Mesh.h"
#include "LightSource.h"
#include "Material.h"

class Scene {
public:
	inline Scene () {}
	virtual ~Scene() {}

	inline const Camera& camera() const { return m_camera; }

	inline Camera& camera() { return m_camera; }
    
    inline const LightSource& lightsource() const { return m_lightsource; }
    
    inline LightSource& lightsource() { return m_lightsource; }
    
    inline const Material& material() const { return m_material; }
    
    inline Material& material() { return m_material; }
	
	inline const std::vector<Mesh> & meshes () const { return m_meshes;  }
	
	inline std::vector<Mesh> & meshes () { return m_meshes; }

private:
	Camera m_camera;
    LightSource m_lightsource;
    Material m_material;
	std::vector<Mesh> m_meshes;
    
};
