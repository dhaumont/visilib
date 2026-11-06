 Result: [Visible, Occluded, Computing, Unknown, Error]

[x] Animate 
        
- Statistics - 
[            ]
[            ]
- ---------- - 

> Light

   Size         0.1 -----+---- 1.0 * 10e-9
   Vertex      [*3,16] +-     
    
   > Gysmo
      Radius        1  ----+----- 10   
      Theta       -90  -----+---- +90
      Phi        -180  -----+---- 180
   
   > Oritentation
      Rotate      -180 ----------180
      Tilt 1      -90  ----------90
      Tilt 2      -90 ----------90
    
> Viewcell

   Size         0.1 -----+---- 1.0 * 10e-9 
   [x] Early stop ON

> Occluders

 - Occluders [Sphere, Cube, slot, Bunny, Sponge, triangle...]      
 - Size       0.1  -----+---- 2.0 * 10e-9
 - Altitude     0   -----+---- 1.0 * 10e-9
    > Geometry
      - Subdivision [1,4] +-
      - Randomness  0  +--------- 50%
      - Holes       0  +--------- 20%
    > Positioning
      - Theta    -90   -----+---- +90
      - Phi     -180   -----+---- 180
      - Offset y  -1.0  -----+---- 1.0 * 10e-9
      - Offset z  -1.0  -----+---- 1.0 * 10e-9
 
 - Ray-tracing: [*Visilib, Embree]

> Receiver

 - Size  0.1 -----+---- 1.0 * 10e-9 
 

> Animation
   [x] Light - Gysmo phi
   [x] Light - Gysmo theta  
   [x] Occluder - Auto-rotate phi
   [x] Occluder - Auto-rotate theta
   [RESET] < > 

> Solver 

[x] Early stop ON

- Epsilon visibility: 1e-6  +-


- Solver: [*Sequential, Aperture, Monte-Carlo]

[Monte Carlo]
   - Number samples: 1024k 

[Sequential, Aperture]
   - [x] Silhouette optimization
   - [x] Local frame split
   - Plucker Normalization [None, 6D, 3D]

> Arithmetic
    - Arithmetic: [Mfpr, GMP rational, double, float]
    - Epsilon arithmetic: 1e-17..1e-15 (hardcoded)
    - Global scaling: 10e[-9,+9] 


> Statistics
[x] Ray splits
[x] Occluders
etc

> Visualization

Results
   [x] Penumbra vertices
   [x] Penumbra edges
   [x] Extremal stabbing lines (free)
   [x] Sampling lines (free)

Geometry
   [x] Light   
   [x] Receiver

Debugging
   [x] Processed edges [*splitted, touched, internal, silhouettes]
   [x] Removed occluders
   [ ] sampling lines (blocked)
   [ ] sampling lines (free)
   [ ] Extremal Stabbing lines (blocked)
   [ ] Extreemal Stabbing lines (degenerated)
   

SAVE   LOAD


-1         0           1    
<----------|----------->   LIGHT  

Altitude




       Occluder
       *******
       *     *
       *     *
       *******






-1         0           1    
<----------|----------->   OCCLUDER