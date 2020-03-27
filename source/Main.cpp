#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <exception>

#include "CommandLine.h"
#include "Image.h"
#include "Ray.h"
#include "Camera.h"
#include "Mesh.h"
#include "Scene.h"
#include "RayTracer.h"
#include "Material.h"
#include "LightSource.h"
#include <chrono>


using namespace std;

void createPlane(Mesh &mesh, vector<Vec3f> corners, Vec3f normal){
    int numVerts = mesh.vertexPositions().size();
    mesh.vertexPositions().insert( mesh.vertexPositions().end(),
                               {corners[0],
                                corners[1],
                                corners[2],
                                corners[3]});
    mesh.vertexNormals().insert( mesh.vertexNormals().end(),
                                {normal,
                                 normal,
                                 normal,
                                 normal});
    mesh.indexedTriangles().push_back(Vec3i(numVerts,numVerts+1, numVerts+3));
    mesh.indexedTriangles().push_back(Vec3i(numVerts,numVerts+2, numVerts+3));
}

int main (int argc, char ** argv) {
	CommandLine args;
	if (argc > 1) {
		try {
			args.parse(argc, argv);
		} catch (const std::exception & e) {
			std::cerr << e.what() << std::endl;
			args.printUsage (argv[0]);
			exit(1);
		}
	}
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
	// Initialization
    
	Image image (args.width (), args.height ());
	Scene scene;
	
	Camera camera(Vec3f(0.3f, 0.2f, 2.3f),
				  Vec3f(),
				  Vec3f(0.f, 1.f, 0.f),
				  60.f,
				  float(args.width()) / args.height());

	scene.camera() = camera;
	
    //LightSource lightsource(Vec3f(0.f, 2.f, 2.f),
    //                        Vec3f(1.f, 1.f, 1.f),
    //                        0.9f);
    
    LightSource lightsource_square(Vec3f(0.3f, 1.f, 1.1f),  // position
                            Vec3f(1.f, 1.f, 1.f),  // color
                            Vec3f(0.f, 0.3f, 0.f),  // direction
                            5.f,                   // intensity
                            1.f),                  // sideLength
    
    lightsource_point(Vec3f(1.3f, 1.f, 0.f),  // position
                        Vec3f(1.f, 1.f, 1.f),  // color
                        Vec3f(0.f, 0.f, 0.f),  // direction
                        5.f,                   // intensity
                        0.001f),                 // sideLength
    
    lightsource_point2(Vec3f(-1.3f, 2.f, 1.f),  // position
                        Vec3f(1.f, 1.f, 1.f),  // color
                        Vec3f(0.f, 0.f, 0.f),  // direction
                        4.7f,                   // intensity
                        0.01f);                 // sideLength
    
    
    
    // define lightsources
    scene.lightsources().push_back (lightsource_square);
    scene.lightsources().push_back (lightsource_point);
    scene.lightsources().push_back (lightsource_point2);
    
    // define materials
    float cube_kd = 0.1f, // diffusion coeff,
          walls_kd = 0.6f,
          head_kd = 0.1f,
          // alpha coeff
          cube_alpha = 0.1f,
          walls_alpha = 0.3f,
          head_alpha = 0.1f;
    // albedo
    Vec3f cube_albedo(0.9f, 0.9f, 0.9f),
          walls_albedo(0.9f, 0.5f, 0.5f),
          head_albedo(0.9f, 0.9f, 0.9f);
    // Fresnel term
    Vec3f head_F0(0.31f, 0.31f, 0.31f), // glass
          cube_F0(1.0f, 0.86f, 0.57f), // gold
          walls_F0(0.5f,0.5f,0.5f);
    
	Material material_head(head_kd, head_alpha, head_albedo, head_F0),
             material_cube(cube_kd, cube_alpha, cube_albedo, cube_F0),
             material_walls(walls_kd, walls_alpha, walls_albedo, walls_F0),
             material_cube2(walls_kd, walls_alpha, Vec3f(0.4f, 0.4f, 0.9f), Vec3f(0.3,0.3,0.3));
    
	Mesh mesh_head, mesh_walls, mesh_cube, mesh_cube2;
    mesh_head.material() = material_head;
    mesh_walls.material() = material_walls;
    mesh_cube.material() = material_cube;
    mesh_cube2.material() = material_cube2;
    //material_cube2;
    
    // Loading meshes
	try {
		mesh_head.loadOFF("../meshes/example_cube.off");
        mesh_cube.loadOFF("../meshes/cube_tri.off");
        mesh_cube2.loadOFF("../meshes/cube_tri2.off");
	}
	catch (const std::exception & e) {
		std::cerr << e.what() << std::endl;
		exit(1);
	}
    
    float box_size = 1.5f;
    // Creating the ground
    vector<Vec3f> corners = {
        Vec3f(box_size,-1.f,box_size),
        Vec3f(box_size,-1.f,-box_size),
        Vec3f(-box_size,-1.f,box_size),
        Vec3f(-box_size,-1.f,-box_size)
    };
    
    Vec3f normal(0.f,1.f,0.f);
    createPlane(mesh_walls,  corners,  normal);
    
    // Creating the back wall
    corners = {
        Vec3f(-box_size,-1.f,-box_size),
        Vec3f(box_size,-1.f,-box_size),
        Vec3f(-box_size,1.25f,-box_size),
        Vec3f(box_size,1.25f,-box_size)
    };
    
    Vec3f normal2(0.f,0.f,1.f);
    createPlane(mesh_walls,  corners,  normal2);
    
    
    mesh_walls.computeBVH();
    mesh_head.computeBVH();
    mesh_cube.computeBVH();
    
    
    scene.meshes ().push_back (mesh_walls);
    scene.meshes ().push_back (mesh_head);
    scene.meshes ().push_back (mesh_cube);
    scene.meshes ().push_back (mesh_cube2);
    
    std::cout << "RayTracer creation: starts";
	RayTracer rayTracer(args.numRays (), args.mode ());
    std::cout << ".....ends." << std::endl;
    
	// Rendering
	std::cout << "Background filling: starts";
	image.fillBackground ();
    std::cout << ".....ends." << std::endl;

    //std::cout << "Ray tracing: starts";
	rayTracer.render (scene, image);
	//std::cout << "ends." << std::endl;
    
    
    
	image.savePPM (args.outputFilename ());
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    std::cout << "Total time is " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
    
    
	

	return 0;
}
