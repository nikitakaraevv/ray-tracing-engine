#pragma once

#include <string>
#include <exception>

/// Helper class to parse the command line arguments
class CommandLine {
public:
	inline CommandLine () :
		m_width(480),
		m_height(270),
        m_numRays(6),
        m_mode(0),
		m_outputFilename("output.ppm") {}
	virtual ~CommandLine() {}

	inline size_t width() const { return m_width; }

	inline size_t height() const { return m_height; }
    
    inline size_t numRays() const { return m_numRays; }
    
    inline size_t mode() const { return m_mode; }

	inline const std::string& outputFilename() const { return m_outputFilename; }

	void printUsage(const char* command) {
		std::cerr << "USAGE: " << command << " [-w/-width <image width>][-h/-height <image height>][-o/-output <outputfilename>][-N/-n/-numRays <number of rays per pixel>][-m/-mode <mode (0 for Ray tracing, 1 for Path tracing, 2 for Path guiding)>]" << std::endl;
	}

	inline void parse(int argc, char** argv) {
		for (int i = 1; i < argc; i++) {
			std::string argi = argv[i];
			if (i == argc - 1) {
				if (argi == "-help") {
					printUsage(argv[0]);
					exit(0);
				}
				else {
					throw std::runtime_error("Missing argument");
				}
			}
			if (argi == "-w" || argi == "-width") {
				m_width = std::atoi(argv[++i]);
			}
			else if (argi == "-h" || argi == "-height") {
				m_height = std::atoi(argv[++i]);
			}
			else if (argi == "-o" || argi == "-output") {
				m_outputFilename = std::string(argv[++i]);
			}
            else if (argi == "-N" || argi == "-n" || argi == "-numRays" ) {
                m_numRays = std::atoi(argv[++i]);
            }
            else if (argi == "-m" || argi == "-mode") {
                m_mode = std::atoi(argv[++i]);
            }
			else {
				throw std::runtime_error(std::string("Unknown argument <" + std::string(argi) + ">").c_str());
			}
		}
	}

private:
	size_t m_width,
    m_height,
    m_numRays,
    m_mode;
	std::string m_outputFilename;
};
