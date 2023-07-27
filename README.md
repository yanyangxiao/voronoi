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

## license FAQ

1) What are the licensing conditions for this project ?

This project is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License (AGPLv3), or (at your option) any later version. A copy of the AGPLv3 is reproduced here. This project is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

2) Sorry, but under these conditions I will be unable to use this project. My company is selling the software it produces and makes a living from that. If we take the AGPL seriously, we would starve, or we cannot use this project.

If the terms and conditions of the AGPL v.3. would prevent you from using this project, please consider the option to obtain a commercial license for a fee. These licenses are provided "as-is", unlimited in time for a one time fee. Please send corresponding requests to yanyangxiaoxyy@gmail.com. Please do not forget to include some description of your company and the realm of its activities.

3) How much is a commercial license then ?

Price qoutes are provided upon request. Please do not hesitate to contact yanyangxiaoxyy@gmail.com for a qoute.
