#pragma once

#include "Image.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

#include "Vec3.h"

inline void Image::fillBackground(const Vec3f& color) {
  for (size_t y = 0; y < m_height; y++)
    for (size_t x = 0; x < m_width; x++) {
      static const Vec3f color0(0.1f, 0.2f, 0.8f);
      static const Vec3f color1(0.9f, 0.9f, 1.0f);
      float alpha =
          std::clamp(static_cast<float>(y) / (m_height - 1), 0.f, 1.f);
      m_pixels[y * m_width + x] = mix(color0, color1, alpha);
    }
}

inline void Image::savePPM(const std::string& filename) {
  std::ofstream out(filename.c_str());
  if (!out) {
    std::cerr << "Cannot open file " << filename.c_str() << std::endl;
    std::exit(1);
  }
  out << "P3" << std::endl
      << m_width << " " << m_height << std::endl
      << "255" << std::endl;
  for (size_t y = 0; y < m_height; y++)
    for (size_t x = 0; x < m_width; x++) {
      out << static_cast<unsigned int>(255.f * m_pixels[y * m_width + x][0])
          << " "
          << static_cast<unsigned int>(255.f * m_pixels[y * m_width + x][1])
          << " "
          << static_cast<unsigned int>(255.f * m_pixels[y * m_width + x][2])
          << " ";
    }
  out << std::endl;
  out.close();
}
