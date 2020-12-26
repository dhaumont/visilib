# Visilib #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

The visilib project aims to provide an efficient and robust C++ library for exact visibility computation in 3D.

The main feature of the current version is exact occlusion query between two polygons, a polygon and a segment, or between two segments, and allows to answer the fundamental question "Are these two polygons mutually visible?". The principal application is from-region occlusion culling and exact Potentially Visible Set computation (PVS). The visibility determination is exact and make no simplification, in opposition to conservative or aggressive visibility methods: visilib performs the complete penumbra fusion of the occluders, by the mean of CSG operations in Plücker space. All the technical details are described in the companion paper "A Low Dimensional Framework for Exact Polygon-to-Polygon Occlusion Queries", presented at EGSR 2005 (cf. download section)

###Features###

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
- fully documented (doxygen)

- Visibility performed via CSG (Constructive Solid Geometry) operations in Plucker space:
     - Optimal minimum polytope building algorithm in Plucker space
     - Polytope splitting algorithm using only skeleton of polytopes (vertices and edges)
     - Silhouette optimization algortihm reducing drastically the number of CSG operations
     - Guided aperture sampling and early termination: rays are casted in the visibility apertures left by the previously processed occluders, leading to early termination in case of mutual visibility
     - Occluder selection using previous queries
     - Computational Geometry predicates in Plucker space

###Applications###
- Potentially Visible Set computation (PVS)
- Sound propagation (edge diffraction)
- Shadow and lighting computation
- Robotic

###Cross-platform###

- CMake: linux, windows

###Supported Compilers ###

- GCC // Add version
- Visual Studio 2019

###Dependencies###
- CMake 3.14.0

###Optional dependencies:###
- CGAL and LEDA (exact arithmetic)
- Embree (fast ray-triangle intersection)
- Google Test (unit testing)
- Doxygen (documentation generation)

###Unit tests###

- based on Google test

###Visilib Extensions###

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
 
