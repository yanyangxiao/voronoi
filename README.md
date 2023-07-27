# voronoi-power-diagram

This project implements the computations of Voronoi diagrams and power diagrams using the Delaunay/regular triangulation based method and the k-nearest neighbors (knn) based method, respectively, while the knn-based method for the computation of power diagrams is proposed in our paper "Yanyang Xiao, Juan Cao, Shaoping Xu and Zhonggui Chen. Meshless Power Diagrams. Computers & Graphics, 2023, (Proc. SMI), 114: 247-256"

## dependencies

CGAL: https://www.cgal.org/ </br>
nanoflann: https://github.com/jlblancoc/nanoflann.git </br>
OpenMesh: http://www.openmesh.org/ (only for loading mesh in the surface case)

## compile

Using CMake

(The include and lib path of OpenMesh in the CMakeLists needs to be adjusted manually to fit your system)

## usage

Once build the code successfully, go to the example/xxx/run-powerx-demo.bat and make sure the path of powerx-demo.exe is correct, then double click the .bat file to run the program, you will obtain 4 result files (powerx-rt.xxx/powerx-rt-omp.xxx/powerx-knn.xxx/powerx-knn-omp.xxx) for each test.

Try to open .svg files or .off files to view the power diagrams.

## cite

If you find our code or paper useful, please consider citing

```
@article{Xiao2023-knnpower,
author = {Xiao, Yanyang and Cao, Juan and Xu, Shaoping and Chen, Zhonggui},
title = {Meshless Power Diagrams},
journal = {Computers & Graphics},
volume = {114},
pages= {247-256},
year = {2023}
}
```