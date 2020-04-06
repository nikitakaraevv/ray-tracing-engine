#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

#include "Vec3.h"

class Image {
 public:
  inline Image(size_t width = 64, size_t height = 64)
      : m_width(width), m_height(height) {
    m_pixels.resize(width * height);
  }

  inline virtual ~Image() {}

  inline size_t width() const { return m_width; }

  inline size_t height() const { return m_height; }

  inline const Vec3f& operator()(size_t x, size_t y) const {
    return m_pixels[y * m_width + x];
  }

  inline Vec3f& operator()(size_t x, size_t y) {
    return m_pixels[y * m_width + x];
  }

  inline void fillBackground(const Vec3f& color = Vec3f(0.f, 0.f, 1.f));

  inline void savePPM(const std::string& filename);

 private:
  size_t m_width;
  size_t m_height;
  std::vector<Vec3f> m_pixels;
};
