#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <exception>

#include "Material.h"
#include "Vec3.h"

typedef Vec3i Triangle;
// A simple container for mesh data
class Mesh {
public:
	const std::vector<Vec3f>& vertexPositions() const { return m_vertexPositions; }

	std::vector<Vec3f>& vertexPositions() { return m_vertexPositions; }

	const std::vector<Vec3f>& vertexNormals() const { return m_vertexNormals; }

	std::vector<Vec3f>& vertexNormals() { return m_vertexNormals; }

	const std::vector<Triangle>& indexedTriangles() const { return m_indexedTriangles; }

	std::vector<Triangle>& indexedTriangles() { return m_indexedTriangles; }
    
    inline const Material& material() const { return m_material; }
    
    inline Material& material() { return m_material; }

	void recomputeNormals() {
		auto& N = vertexNormals();
		const auto& P = vertexPositions();
		const auto& T = indexedTriangles();
		N.resize(P.size(), Vec3f());
		for (size_t i = 0; i < T.size(); i++) {
			Vec3f nt(triangleNormal(i));
			for (size_t j = 0; j < 3; j++)
				N[T[i][j]] += nt;
		}
		for (size_t i = 0; i < N.size(); i++)
			N[i].normalize();
	}

	inline void loadOFF(const std::string& filename) {
		try {
			m_vertexPositions.clear();
			m_indexedTriangles.clear();
			std::ifstream in(filename.c_str());
			if (!in)
				throw std::runtime_error(("Error loading OFF file: " + filename).c_str());
			std::string offString;
			unsigned int numV, numF, numE;
			in >> offString;
			skipHashCommentLine(in);
			in >> numV >> numF >> numE;
			skipHashCommentLine(in);
			auto& P = vertexPositions();
			P.resize(numV);
			for (unsigned int i = 0; i < numV; i++)
				in >> P[i];
			auto& T = indexedTriangles();
			unsigned int s;
			for (unsigned int i = 0; i < numF; i++) {
				in >> s;
				std::vector<unsigned int> v(s);
				for (unsigned int j = 0; j < s; j++)
					in >> v[j];
				for (unsigned int j = 2; j < s; j++)
					T.push_back(Triangle(v[0], v[j - 1], v[j]));
			}
			in.close();
		}
		catch (const std::exception & e) {
			throw std::runtime_error((std::string("Error Loading OFF file: ") + std::string(e.what())).c_str());
		}
		recomputeNormals();
	}

private:
	Vec3f triangleNormal(size_t i) const {
		return normalize(cross(m_vertexPositions[m_indexedTriangles[i][1]] - m_vertexPositions[m_indexedTriangles[i][0]],
							   m_vertexPositions[m_indexedTriangles[i][2]] - m_vertexPositions[m_indexedTriangles[i][0]]));
	}

	static const unsigned int MAX_COMMENT_SIZE = 1024;

	void skipHashCommentLine(std::ifstream& in) {
		while (in.peek() == '\n' || in.peek() == ' ') // remove former separators
			in.get();
		int c = in.peek();
		if (c == '#') {
			char trash[MAX_COMMENT_SIZE];
			in.getline(trash, MAX_COMMENT_SIZE - 1);
		}
	}

	std::vector<Vec3f> m_vertexPositions;
	std::vector<Vec3f> m_vertexNormals;
	std::vector<Triangle> m_indexedTriangles;
    Material m_material;
};
