#define _USE_MATH_DEFINES
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "Camera.h"
#include "CommandLine.h"
#include "Image.cpp"
#include "LightSource.h"
#include "Material.h"
#include "Mesh.h"
#include "Ray.h"
#include "RayTracer.h"
#include "Renderer.h"
#include "Scene.h"

using namespace std;

inline void createPlane(Mesh &mesh, vector<Vec3f> corners, Vec3f normal) {
  int numVerts = mesh.vertexPositions().size();
  mesh.vertexPositions().insert(
      mesh.vertexPositions().end(),
      {corners[0], corners[1], corners[2], corners[3]});
  mesh.vertexNormals().insert(mesh.vertexNormals().end(),
                              {normal, normal, normal, normal});
  mesh.indexedTriangles().push_back(
      Vec3i(numVerts, numVerts + 1, numVerts + 3));
  mesh.indexedTriangles().push_back(
      Vec3i(numVerts, numVerts + 2, numVerts + 3));
}

void createCornellBox(float box_size, float ceiling, Mesh &mesh_walls,
                      Mesh &mesh_right_wall, Mesh &mesh_left_wall) {
  // Creating the ground
  vector<Vec3f> corners = {
      Vec3f(box_size, -1.f, box_size), Vec3f(box_size, -1.f, -box_size),
      Vec3f(-box_size, -1.f, box_size), Vec3f(-box_size, -1.f, -box_size)};

  Vec3f normal(0.f, 1.f, 0.f);
  createPlane(mesh_walls, corners, normal);

  // Creating the back wall
  corners = {Vec3f(-box_size, -1.f, -box_size),
             Vec3f(box_size, -1.f, -box_size),
             Vec3f(-box_size, ceiling, -box_size),
             Vec3f(box_size, ceiling, -box_size)};

  normal = {0.f, 0.f, 1.f};
  createPlane(mesh_walls, corners, normal);

  // Creating the ceiling

  corners = {Vec3f(box_size, ceiling, box_size),
             Vec3f(box_size, ceiling, -box_size),
             Vec3f(-box_size, ceiling, box_size),
             Vec3f(-box_size, ceiling, -box_size)};

  normal = {0.f, -1.f, 0.f};
  createPlane(mesh_walls, corners, normal);

  // Creating the left wall

  corners = {Vec3f(-box_size, -1.f, box_size),
             Vec3f(-box_size, -1.f, -box_size),
             Vec3f(-box_size, ceiling, box_size),
             Vec3f(-box_size, ceiling, -box_size)};

  normal = {1.f, 0.f, 0.f};
  createPlane(mesh_left_wall, corners, normal);

  // Creating the right wall

  corners = {Vec3f(box_size, -1.f, box_size), Vec3f(box_size, -1.f, -box_size),
             Vec3f(box_size, ceiling, box_size),
             Vec3f(box_size, ceiling, -box_size)};

  normal = {-1.f, 0.f, 0.f};
  createPlane(mesh_right_wall, corners, normal);
}

inline void rotationY(Mesh &mesh, float phi) {
  Vec3<Vec3f> transformMatrix(Vec3f(cos(phi), 0, sin(phi)), Vec3f(0, 1, 0),
                              Vec3f(-sin(phi), 0, cos(phi)));
  vector<Vec3f> vertexPositions(mesh.vertexPositions());
  for (Vec3f &position : vertexPositions) {
    Vec3f newPosition(dot(transformMatrix[0], position),
                      dot(transformMatrix[1], position),
                      dot(transformMatrix[2], position));
    position = newPosition;
  }
  mesh.vertexPositions() = vertexPositions;
}

void initLightSources(Scene &scene) {
    // Initialize lightsources
    LightSource lightsource_point(Vec3f(-1.4f, 1.f, 2.9f),  // position
                                  Vec3f(1.f, 1.f, 1.f),     // color
                                  Vec3f(0.3f, 0.f, -1.f),   // direction
                                  0.85f,                      // intensity
                                  0.01f),                   // sideLength

        lightsource_point2(Vec3f(1.4f, 1.f, 2.9f),   // position
                           Vec3f(1.f, 1.f, 1.f),     // color
                           Vec3f(-0.3f, 0.f, -1.f),  // direction
                           0.85f,                      // intensity
                           0.01f),                   // sideLength

        lightsource_square(Vec3f(0.f, -0.3f, 1.1f),  // position
                           Vec3f(1.f, 1.f, 1.f),     // color
                           Vec3f(0.f, 0.f, -1.f),    // direction
                           0.85f,                      // intensity
                           0.1f);                    // sideLength

    scene.lightsources().push_back(lightsource_point);
    scene.lightsources().push_back(lightsource_point2);
    scene.lightsources().push_back(lightsource_square);

}

void initLightMaterials(Mesh & mesh_walls,
                        Mesh & mesh_cube,
                        Mesh & mesh_cube2,
                        Mesh &mesh_left_wall,
                        Mesh & mesh_right_wall) {
    // Initialize materials of the meshes
    float cube_kd = 0.1f,  // diffusion coeff,
        walls_kd = 0.6f,
          // alpha coeff
        cube_alpha = 0.1f, walls_alpha = 0.3f;
    // albedo
    Vec3f cube_albedo(0.9f, 0.9f, 0.9f), walls_albedo(0.96f, 0.96f, 0.86f);
    // Fresnel term
    Vec3f //(0.31f, 0.31f, 0.31f)  // glass
        cube_F0(1.0f, 0.86f, 0.57f),     // gold
        walls_F0(0.5f, 0.5f, 0.5f);

    Material material_cube(cube_kd, cube_alpha, cube_albedo, cube_F0),
        material_walls(walls_kd, walls_alpha, walls_albedo, walls_F0),
        material_lw(walls_kd, walls_alpha, Vec3f(0.9f, 0.3f, 0.3f), walls_F0),
        material_rw(walls_kd, walls_alpha, Vec3f(0.3f, 0.9f, 0.3f), walls_F0),
        material_cube2(0.8f, 0.9f, Vec3f(0.4f, 0.4f, 0.9f), Vec3f(0.3, 0.3, 0.3));

    
    mesh_walls.material() = material_walls;
    mesh_cube.material() = material_cube;
    mesh_cube2.material() = material_cube2;
    mesh_left_wall.material() = material_lw;
    mesh_right_wall.material() = material_rw;
}

int main(int argc, char **argv) {
  CommandLine args;
  if (argc > 1) {
    try {
      args.parse(argc, argv);
    } catch (const std::exception &e) {
      std::cerr << e.what() << std::endl;
      args.printUsage(argv[0]);
      exit(1);
    }
  }
  chrono::steady_clock::time_point begin = chrono::steady_clock::now();
  // Initialization
  Image image(args.width(), args.height());
  Scene scene;

    
  // Initialize camera
  Camera camera(Vec3f(0.3f, 0.6f, 2.3f), Vec3f(), Vec3f(0.f, 1.f, 0.f), 60.f,
                float(args.width()) / args.height());

  scene.camera() = camera;
    
  // Initialize lightsources
  initLightSources(scene);
    
  // Initialize meshes
  Mesh mesh_walls, mesh_cube, mesh_cube2, mesh_left_wall,
        mesh_right_wall;
    
  // Initialize materials of the meshes
  initLightMaterials( mesh_walls,
                      mesh_cube,
                      mesh_cube2,
                      mesh_left_wall,
                      mesh_right_wall);

  // Loading meshes
  try {
    mesh_cube.loadOFF("../meshes/cube_tri.off");
    mesh_cube2.loadOFF("../meshes/cube_tri2.off");
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    exit(1);
  }

  // Creating the Cornell box scene
  float box_size = 1.51f, ceiling = 1.5f;
  createCornellBox(box_size, ceiling, mesh_walls, mesh_right_wall,
                   mesh_left_wall);
  // mesh_walls.computeBVH();
  // mesh_cube.computeBVH();
    
  // Rotate cubes in the scene
  rotationY(mesh_cube, M_PI / 4.5f);
  rotationY(mesh_cube2, -M_PI / 4.5f);

  scene.meshes().push_back(mesh_walls);
  scene.meshes().push_back(mesh_left_wall);
  scene.meshes().push_back(mesh_right_wall);
  scene.meshes().push_back(mesh_cube);
  scene.meshes().push_back(mesh_cube2);
    
  // Create different renderers depending on usage of photon mapping
  RayTracer rayTracer;
  Renderer renderer;
  if (args.numPhotons() > 0)
    renderer = Renderer(scene, args.numRays(), args.mode(), rayTracer,
                        args.numPhotons(), args.k());
  else
    renderer = Renderer(scene, args.numRays(), args.mode(), rayTracer);

  image.fillBackground();

   // Rendering
  renderer.render(image);
    
  // Save final ppm image
  image.savePPM(args.outputFilename());
  chrono::steady_clock::time_point end = chrono::steady_clock::now();
  std::cout
      << "Total time is "
      << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count()
      << "[s]" << std::endl;

  return 0;
}
