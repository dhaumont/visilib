### Visilib Overview ###

The visilib project aims to provide an efficient and robust library for **exact visibility computation** in 3D.

The main feature of the current version of Visilib is the exact visibility computation between two primitives (points, segments or polygons). The visibility computation is exact and make no simplification, in opposition to conservative or aggressive visibility methods: Visilib performs the complete penumbra fusion of the occluders, by the use of CSG operations in Plücker space.

The algorithm is described in detail in the technical paper "A Low Dimensional Framework for Exact Polygon-to-Polygon Occlusion Queries" (see **Paper**).
However, the library has been completely rewritten and benefits from several improvements in comparison to the paper (see **Features**).

Visilib is written in C++, it is **cross-platform** and distributed as a **header only library**.


### Features ###

- Available exact visibility query:
    - polygon to polygon
    - segment to segment
    - segment to polygon
    - point to point
    - point to segment
    - point to polygon
    
- robust arithmetic computation via exact artithmetic via CGAL and LEDA libraries  (multiple precision arithmetic combined with interval arithmetic computations)
- fast ray-triangle intersection based on Intel Embree library
- easy to deploy: header only library 

- Visibility performed via CSG (Constructive Solid Geometry) operations in Plucker space:
     - Optimal minimum polytope building algorithm in Plucker space
     - Polytope splitting algorithm using only skeleton of polytopes (vertices and edges)
     - Silhouette optimization algortihm reducing drastically the number of CSG operations
     - Guided aperture sampling and early termination: rays are casted in the visibility apertures left by the previously processed occluders, leading to early termination in case of mutual visibility
     - Occluder selection using previous queries
     - Computational Geometry predicates in Plucker space

### Applications ###
- Potentially Visible Set computation (PVS)
- Sound propagation (edge diffraction)
- Shadow and lighting computation
- Robotic

### Cross-platform ###

- CMake: linux, windows

### Supported Compilers ###

- GCC // Add version
- Visual Studio 2019

### Dependencies ###
- CMake 3.14.0

### Optional dependencies ###
- CGAL and LEDA (exact arithmetic)
- Embree (fast ray-triangle intersection)
- Google Test (unit testing)
- Doxygen (documentation generation)

## #Unit tests ###

- based on Google test

### Visilib Extensions ###

Visilib will be extended in the future, in order to provide a complete framework for exact visibility determination in 3D. The following features will be added in the future:

 TODO 
 Full visibility determination (in addition to the boolean visible/hidden predicates) for exact form-factor computation and exact shadow computation

Project Information

    License: GNU GPL v3

Labels:
Algorithm 3D graphics geometry Mathematics Plucker Coordinates Library pvs occlusion visibility exact culling computation CPlusPlus CSG


* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up
TODO 
* Configuration
TODO 
* Dependencies
TODO 
* Database configuration
TODO 
* How to run tests
TODO 
* Deployment instructions
TODO 

### Contribution guidelines ###

* Writing tests
TODO 
* Code review
TODO 
* Other guidelines
TODO 

### Who do I talk to? ###

* Haumont Denis (denis.haumont@gmail.com)


### License ###
Visilib is  GNU GPL v3 licensed.
If you use Visilib in a publication, please cite the following paper and refer to this project.
```
@inproceedings{
  title     = {A Low Dimensional Framework for Exact Polygon-to-Polygon Occlusion Queries},
  author    = {Denis Haumont and
               Otso Makinen and
               Shaun Nirenstein}, 
  booktitle = {Proceedings of the Eurographics Symposium on Rendering Techniques,
               Konstanz, Germany, June 29 - July 1, 2005},
  pages     = {211--222},
  publisher = {Eurographics Association},
  year      = {2005}
}
```


### visilibDemo ###

visilibDemo is an interactive demo application provided with visilib library.

It allows to test the different library options interactively, by pressing the **following keys**:

- 's': enable/disable silhouette optimisation

Without silhouette optimisation, all the edges of the occluders are used to perform CSG operations in Plucker space.
With silhouette optimisation, only the potential silhouette edges of the occluders with respect to the source polygons are used.

- 'n': enable/disable nomalization

When normalisation is activated, all points in the projective Plucker space are projected onto an unit hypersphere, in order to avoid numerical errors.

- 'r': enable/disable representative line sampling strategy
  
  
- 'f': enable/disable detect aperture only

When this option is activated, the visibility computation stops as soon as the two sources are proven to be mutually visible. Otherwise, the computation continues until the complete visibility information has been computed.
 
 
- 'f': enable/disable fast silhouette rejection test
 
- '+/-': increase/decrease scaling of query polygons



- '1/2': increase/decrease number of vertices of the first query source

The  key '1' decreases the number of vertices of the first source.
The  key '2' increases the number of vertices of the first source.

- 'w': write config

- 'o': open config

- 'space': start/pause animation

When animation is started, the source polygons will perform a rotation arount the scene.

- 'Enter': show/hide geometry

For debugging purpose, allow to change the representation of the scene and of the query polygons. 

'h': display help
Display the list of keys and their effect


WIndows (VS2019 x64)


1) EXACT ARITHMETIC: CGAL + BOOST + LEDA

CGAL-4.13.1
   1. download and install BOOST version 1_70_0
https://www.boost.org/
   
   2. define BOOST_INCLUDEDIR variable
   
Note: you can decomment the line in generateProject.bat
   
   3. download and install CGAL version 4.13.1
   https://github.com/CGAL/cgal/releases
   4. define CGAL_DIR envionment variable
Note: you can decomment the line in generateProject.bat
 

	5. download and install LEDA version 6.3
	http://www.algorithmic-solutions.com/index.php/products/leda-free-edition
 
 
   5. activate   


2) FAST RAY TRACING: EMBREE


Embree 3.0
   1. download and install embree
     https://www.embree.org/downloads.html
   2. define EMBREE_INSTALL_DIR environment variable
       in the compilation command, call 
	   set EMBREE_INSTALL_DIR=c:\Program Files\Intel\Embree3 x64
   
   Note: you can decomment the line in generateProject.bat
   
   3. activate "embree" flag in CMake
   
   Open "CMakeList.txt"
   set (USE_EMBREE True)

   4. run generateProject.bat
   
   
3) ADDITIONAL DLL FOR THE DEMO
 - freeGlut
 - Opengl
 
