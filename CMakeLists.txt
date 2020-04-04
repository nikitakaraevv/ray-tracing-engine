cmake_minimum_required(VERSION 3.9)

project(MyRayTracer LANGUAGES CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)


file(GLOB SRCFILES *.cpp)
add_executable (
	MyRayTracer
	source/Main.cpp
	source/CommandLine.h
	source/Vec3.h
	source/Ray.h
	source/Ray.cpp
	source/Image.h
	source/Camera.h
	source/Mesh.h
	source/Scene.h
	source/RayTracer.h
	source/Material.h
	source/LightSource.h
)

set_target_properties(MyRayTracer PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)



