## Visilib Overview 

The Visilib project aims to provide an efficient and robust library for **exact visibility computation** in 3D.

Given a polygonal scene as input treated as occluders, the library allows to determine if two primitives called sources are mutually visible. In case of mutual visibility, Visilib optionally computes the set of lines stabbing the two primitives, providing a complete description of the mutual visibility information. 

The visibility computations are computed in an **exact** way, without making any simplification.

The core algorithm is based on CSG operations in a 5D projective space called Pücker space. If you are interested in implementation details please refer to this technical [paper](https://github.com/dhaumont/visilib/blob/main/paper/paper_egsr_2005.pdf) describing the core algorithm.
However, the library has been completely rewritten and benefits from several improvements (see [Features](https://github.com/dhaumont/visilib#Features)).

Visilib is written in C++, it is **cross-platform** and distributed as a **header only library**.

## Features 

- Visibility computation is exact:  contrary to conservative or aggressive visibility methods, the Visilib library does not make any assumption and performs the complete penumbra fusion of the occluders 
- Visibility performed via CSG (Constructive Solid Geometry) operations in Plücker space:
     - Polytope splitting algorithm using only skeleton of polytopes (vertices and edges)
     - Silhouette optimization algortihm reducing drastically the number of CSG operations
     - Guided aperture sampling and early termination: rays are casted in the visibility apertures left by the previously processed occluders, leading to early termination in case of mutual visibility
     - Occluder selection using previous queries
     - Computational Geometry predicates in Plucker space

**NEW (since paper publication)**

- The source primitives that can be used are polygons, segments and points, allowing to perform polygon-to-polygon, polygon-to-segment, polygon-to-point and point-to-segment visibility query
- Optimal minimum polytope building algorithm in Plucker space (TODO refer to paper),
- Robust arithmetic computation handling based CGAL and LEDA libraries  (multiple precision arithmetic combined with interval arithmetic computations)
- fast ray-triangle intersection based on Intel Embree library
- easy to deploy: header only library

## Applications 
- Potentially Visible Set computation (PVS)
- Sound propagation (edge diffraction)
- Shadow and lighting computation

## Dependencies 
- [CMake] (https://cmake.org/) for cross-platform compilation

## Optional dependencies 
- [CGAL] (https://www.cgal.org/) and [Leda] for exact arithmetic support
- [Intel Embree] (https://www.embree.org), for fast ray-tracing
- [Doxygen] (https://www.doxygen.nl/index.html) for documentation generation
- [freeglut] (http://freeglut.sourceforge.net/) and [OpenGl] (https://www.opengl.org/) for the interactive demo application

## #Unit tests 

## Visilib Extensions 

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

## How do I get set up? 

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



## Who do I talk to? 

* Haumont Denis (denis.haumont@gmail.com)


## License 
Visilib is [GNU GPL v3](https://www.gnu.org/licenses/gpl-3.0) licensed
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


## visilibDemo 

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
 
