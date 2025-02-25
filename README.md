# Hilbert Curve Projection for Optimal Transport Sampling on a Sphere Surface


This repository aims to generate blue noise samples on the sphere by Hilbert Curve Projection.


![sphere_noise](out/merged_sphere_noise.jpg)


(regular, whitenoise, dartthrowing, stratified, poissondisk, NESOTS, HCPSOTS)


![mesh](out/merged_mesh.jpg)
The first row is the result of NESOTS, the second row is the result of HCPSOTS.

<!-- ![](https://github.com/boxaio/HCPSOTS/blob/main/out/merged_sphere_noise.jpg) -->

This repository is based on the following repositories :

[NESOTS](https://github.com/baptiste-genest/NESOTS)

[CEPS](https://github.com/MarkGillespie/CEPS.git)

##  Getting started

We use CMake to build the project with the following commands from the root folder :
 ```bash
 mkdir build 
 cd build
 cmake ..
 make -j20
 ```

## Usage
In the build folder, you can run the following command to generate samples on a sphere :
```bash
./tests
```

## Third Party Libraries
- [geometry-central](https://libigl.github.io/) 
- [polyscope](https://polyscope.run/) 
- [eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) 
- [bvh](https://github.com/madmann91/bvh)


