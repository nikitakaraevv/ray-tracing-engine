cmake_minimum_required(VERSION 3.9)

project(RayTracer LANGUAGES CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)


file(GLOB SRCFILES *.cpp)
add_executable (
	RayTracer
	source/Main.cpp
	source/CommandLine.h
	source/Ray.cpp
)

set_target_properties(RayTracer PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)
