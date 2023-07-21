## Visilib Overview 

The Visilib project aims to provide an efficient and robust library for **exact visibility computation** in 3D.

Given a set of polygonal geometries that are treated as occluders, the library allows to determine if two primitives called sources are mutually visible. In case of mutual visibility, Visilib optionally computes the set of lines stabbing the two primitives, providing a complete description of the mutual visibility information. 

The from-region visibility computations are computed in an **exact** way, without making any simplification.

The core algorithm is based on CSG (Constructive Solid Geometry) operations in a 5D projective space called Pücker space. If you are interested in implementation detail please refer to this technical [paper](https://github.com/dhaumont/visilib/blob/main/paper/paper_egsr_2005.pdf) describing the core algorithm.
However, since this publication the library has been completely rewritten and benefits from several improvements (see [Features](https://github.com/dhaumont/visilib#Features)).

Visilib is written in C++, it is **cross-platform** and distributed as a **header only library**.

## Features 

- From-region visibility computations are exact:  contrary to conservative or aggressive visibility methods, the Visilib library does not make any assumption and performs the complete penumbra fusion of the occluders 
- Visibility performed via CSG (Constructive Solid Geometry) operations in Plücker space:
     - Polytope splitting algorithm using only skeleton of polytopes (vertices and edges)
     - Silhouette optimization algortihm reducing drastically the number of CSG operations
     - Guided aperture sampling and early termination: rays are casted in the visibility apertures left by the previously processed occluders, leading to early termination in case of mutual visibility
     - Occluder selection using already computed visibility information
     - Computational geometry predicates in Plücker space

**NEW (since paper publication)**

- The source primitives that can be used are polygons, segments and points, allowing to perform polygon-to-polygon, polygon-to-segment, polygon-to-point, and point-to-segment visibility query
- Optimal minimum polytope building algorithm in Plucker space (as described in [Charneau et al.]( https://www.cs.purdue.edu/cgvlab/popescu/Visibility/Additional%20References/Charneau%202007.pdf)),
- Robust arithmetic computation based on [CGAL](https://www.cgal.org/) and [Leda](https://www.algorithmic-solutions.com/index.php/products/leda-for-c) libraries  (multiple precision arithmetic combined with interval arithmetic computations)
- Fast ray tracing based on [Intel Embree](https://www.embree.org) library
- Easy to deploy: header-only library
- Filtering of small polytopes (next release)
- Early stop or complete skeleton extraction
- Watertight ray-triangle intersection
- Pre-filtering of useless CGS splits
- Output sensitive occluder selection using optimized convex hull construction
- Lazy on-demand computation of occluders of silhouettes
- Representative line sampling computation correction
   
## Applications
- Potentially Visible Set computation (PVS)
- Sound propagation (edge diffraction)
- Shadow and lighting computation

## Dependencies

The only mandatory dependency is [CMake](https://cmake.org/) used for cross-platform compilation.

The following dependencies are optional: 
- [CGAL](https://www.cgal.org/) and [Leda](https://www.algorithmic-solutions.com/index.php/products/leda-for-c) and [Boost](https://www.boost.org/) for optional exact arithmetic support
- [Intel Embree](https://www.embree.org), for optional fast ray-tracing support
- [freeglut](http://freeglut.sourceforge.net/) and [OpenGl](https://www.opengl.org/) for the interactive demo application

We also use [Doxygen](https://www.doxygen.nl/index.html) for documentation generation.
##  Unit tests 

Visilib is partially unit tested, using the CMake framework. The implementation of the unit test is provided in the [/test/src/](https://github.com/dhaumont/visilib/tree/main/test/src) folder. 
The name of executable to run the unit tests is **visilibTest**.

## How do I get set up? 

* Summary of set up

LINUX 
```
#install required dependency
sudo apt-get install freeglut3-dev
sudo apt-get install cmake

#clone repository
git clone https://github.com/dhaumont/visilib.git

#configure 
cd visilib
source configure.sh

#build
cd build
make -j10

#Launch demo
demo/visilibDemo 
```

WINDOWS


## Visilib License 
Visilib is [GNU GPL v3](https://www.gnu.org/licenses/gpl-3.0) licensed.

If you use Visilib in a publication, please refer to the project homepage ([https://github.com/dhaumont/visilib](https://github.com/dhaumont/visilib)) and  cite the following paper:

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


##  Interactive Demo 

visilibDemo is an interactive demo application provided with visilib library.

It allows to test the different library options interactively, by pressing the **following keys**:

| Key | Usage | Description |
|-----|-------|-------------|
| f | enable/disable detect aperture only| When activated, the visibility computation stops as soon as the two sources are proven to be mutually visible. Otherwise, the computation continues until the complete visibility information has been computed |
| s | enable/disable silhouette optimisation |When activated, only the potential silhouette edges of the occluders wrt the source polygons are considered for CSG operations |
| n | enable/disable hypersphere nomalization| When activated, all the Plücker points are projected onto an unit hypersphere, in order to reduce numerical errors |
| r | enable/disable representative line sampling strategy| |
| f | enable/disable fast silhouette rejection test| |
| + | increase source size | Increase the size of the source polygons |
| - | decrease source size | Decrease the size of the source polygons |
| 1 | decrease source vertex count | Decrease the number of vertices of the first source polygon |
| 2 | increase source vertex count | Increase the number of vertices of the first source polygon |
| w | write config|
| o | open config|
|space| start/pause sources animation | When animation is started, the source polygons will perform a rotation around the scene |
|Enter| show/hide geometry | For debugging purpose, allow to change the representation of the scene and of the query polygons |
|h| display help | Display the list of keys and their effect |


## Acknowledgements

We would like to thanks the following persons for their contribution to this library:
- Otso Mäkinen for the silhouette computation implementation
- Shaun Nirenstein for all the technical discussions
- Eric Lengyel for the polygon clipping algorithm
- Joseph O'Rourke and Cédric Laugerotte for the convex hull implementation
- Sven Woop et al. for the watertight triangle intersection code
- Bajaj and Patucci for the interpolation algorithm code
- David Eberly for the AABB to ray intersection code
- Jiry Bittner for the Plücker point anchor computation
- Meshlab team for the obj reader code

We used several useful libraries and tools to implement Visilib (see [Dependencies](https://github.com/dhaumont/visilib#Dependencies)). We would like to especially thank their authors for their great work and publishing the code.


## Who do I talk to? 

For any question about this library, please contact Denis Haumont.


