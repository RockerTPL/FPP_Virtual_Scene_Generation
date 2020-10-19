# Functioning

This project is a study and realization of the generation of a first-person perspective virtual scene.



# Files

You can find  folders listed below.

- **1_projection_world**
  - **projection_world.m** : Imaging based on the world coordinates.
  - **projection_world_triangle.m** : Imaging based on the world coordinates, where triangles are used to handle the sheltering relationship.
- **2_projection_pixel**
  - **projection_pixel.m** : Imaging based on the pixel plane.
- **3_motion_scene**
  - **motion_edge_projection.m** : Generation of a motion scene where edges are extracted by projection from the camera coordinate system.
  - **motion_edge_projection_c.m** : Generation of a motion scene where edges are extracted by projection from the camera coordinate system, which is realized with C++.
  - **motion_sobel.m** : Generation of a motion scene where edges are extracted by sobel operators. (This result is used in the report and the presentation of the project)
  - **getEg.cpp** : Edge extraction with projection from the camera coordinate system.
  - **bigEg.cpp** : Edge expansion with the distance map.
  - **getNeibTr.cpp** : Find nearby triangles while updating a pixel in the motion scene.
  - **travTr.cpp** : Calculation of the relationship between a pixel and the triangles detected
  - **outil.h** : Basic funtions used in the C++ files.



# Application

- You can simply run the `.m` files in MATLAB to get the results.
- **ATTENTION** : If you use the **Windows** operating system, please use `mex` to compile (in MATLAB) the C++ files in the folder **3_motion_scene** before running the `.m` files in this folder.