# voronoi

This program is the implementation of the paper "Xiao, Yanyang, Cao, Juan, Xu, Shaoping and Chen, Zhonggui. Meshless Power Diagrams. Computers & Graphics, 2023, accept"

## dependencies

CGAL: https://www.cgal.org/ </br>
nanoflann: https://github.com/jlblancoc/nanoflann.git </br>
OpenMesh: http://www.openmesh.org/ (only for loading mesh in the surface case)

## compile

Using CMake

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
volume = {accept},
year = {2023}
}
```