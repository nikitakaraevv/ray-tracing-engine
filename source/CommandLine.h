#pragma once

#include <exception>
#include <string>

/// Helper class to parse the command line arguments
class CommandLine {
 public:
  inline CommandLine()
      : m_width(480),
        m_height(270),
        m_numRays(6),
        m_mode(0),
        m_numPhotons(0),
        m_k(5),
        m_outputFilename("output.ppm") {}
  virtual ~CommandLine() {}

  inline size_t width() const { return m_width; }

  inline size_t height() const { return m_height; }

  inline size_t numRays() const { return m_numRays; }

  inline size_t mode() const { return m_mode; }

  inline size_t numPhotons() const { return m_numPhotons; }

  inline size_t k() const { return m_k; }

  inline const std::string& outputFilename() const { return m_outputFilename; }

  void printUsage(const char* command) {
    std::cerr
        << "USAGE: " << command
        << " [-w/-width <image width>][-h/-height <image height>][-o/-output "
           "<outputfilename>][-N/-n/-numRays <number of rays per "
           "pixel>][-m/-mode <mode (0 for Ray tracing, 1 for Path tracing)>][-p/-numPhotons <number of photons for a photon map. If defined, photon map-based rendering is used.>][-k <number of neighbours in photon mapping. Use only with -p/-numPhotons>]"
        << std::endl;
  }

  inline void parse(int argc, char** argv) {
    for (int i = 1; i < argc; i++) {
      std::string argi = argv[i];
      if (i == argc - 1) {
        if (argi == "-help") {
          printUsage(argv[0]);
          exit(0);
        } else {
          throw std::runtime_error("Missing argument");
        }
      }
      if (argi == "-w" || argi == "-width") {
        m_width = std::atoi(argv[++i]);
      } else if (argi == "-h" || argi == "-height") {
        m_height = std::atoi(argv[++i]);
      } else if (argi == "-o" || argi == "-output") {
        m_outputFilename = std::string(argv[++i]);
      } else if (argi == "-N" || argi == "-n" || argi == "-numRays") {
        m_numRays = std::atoi(argv[++i]);
      } else if (argi == "-m" || argi == "-mode") {
        m_mode = std::atoi(argv[++i]);
      } else if (argi == "-p" || argi == "-numPhotons") {
        m_numPhotons = std::atoi(argv[++i]);
      } else if (argi == "-k") {
        m_k = std::atoi(argv[++i]);
      } else {
        throw std::runtime_error(
            std::string("Unknown argument <" + std::string(argi) + ">")
                .c_str());
      }
    }
  }

 private:
  size_t m_width, m_height, m_numRays, m_mode, m_numPhotons, m_k;
  std::string m_outputFilename;
};
