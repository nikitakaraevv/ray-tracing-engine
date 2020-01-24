#pragma once

#include <string>
#include <exception>

/// Helper class to parse the command line arguments
class CommandLine {
public:
	inline CommandLine () :
		m_width(480),
		m_height(270),
		m_outputFilename("output.ppm") {}
	virtual ~CommandLine() {}

	inline size_t width() const { return m_width; }

	inline size_t height() const { return m_height; }

	inline const std::string& outputFilename() const { return m_outputFilename; }

	void printUsage(const char* command) {
		std::cerr << "USAGE: " << command << " [-w/-width <image width>][-h/-height <image height>][-o/-output <outputfilename>" << std::endl;
	}

	inline void parse(int argc, char** argv) {
		for (int i = 1; i < argc; i++) {
			std::string argi = argv[i];
			if (i == argc - 1) {
				if (argi == "-h" || argi == "-help") {
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
			else {
				throw std::runtime_error(std::string("Unknown argument <" + std::string(argi) + ">").c_str());
			}
		}
	}

private:
	size_t m_width;
	size_t m_height;
	std::string m_outputFilename;
};
