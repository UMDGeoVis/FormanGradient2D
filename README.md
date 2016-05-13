##Forman Gradient 2D
#####A comprehensive library for computing a Forman Gradient on triangle meshes

---
####Library Description
***

The library is composed by two main parts:
- the C++ library collecting all the basic functions for managing the simplicial mesh (LibMesh)
- the C++ library collecting all the basic functions for managing the Forman gradient mesh (LibForman)

For a brief description of each function look at LibMesh/Mesh.h or LibForman/formangradient.h



---
####Compile
***

The project can be compiled with QTCreator or using the terminal by typing “qmake” and “make” from this folder.

The code has been tested under Linux and Mac only. No external libraries are needed.




---
####Furnished main
***

The file main.cpp presents some of the basic functionalities of the library.

1. Read a Mesh file (```Reader::readOFFMesh(mesh,meshFile);```)

2. Read a Scalar field (```Reader::readScalarField(mesh,fieldFile);```)

3. Computer the Forman Gradient
    ```
    FormanGradientVector gradient = FormanGradientVector(&mesh);
    gradient.homotopy_expansion();
    ```

The input parameters expected are:

1. Type of input
  * -i for files provided of scalar field
  * -f if the scalar field is in a separate input file
2. Input files
3. Algorithm for computing the Forman gradient
  * [exp] for using Simulation of Simplicity in combination of the algorithm by Homotopy Expansion
  * [exp-improved] for using he improved Simulation of Simplicity in combination of the algorithm by Homotopy Expansion
  * [immersion] for computing the Forman gradient starting from a watershed decomposition


Examples:

1. Read the terrain data from Crater.off and compute the Forman gradient by homotopy expansion

    ```
    ./watershedToForman -i Crater.off exp
    ```

2. Read the mesh Rim01.off and the scalar field Rim01Field.fld and compute the Forman gradient using the watershed decomposition

    ```
    ./watershedToForman -f Rim01.off Rim01Field.fld immmersion
    ```

---
####Attribution
***

If you use FormanGradient2D in you project please cite the original paper
```
@article{Weiss13,
  author    = {Kenneth Weiss and
               Federico Iuricich and
               Riccardo Fellegara and
               Leila De Floriani},
  title     = {A primal/dual representation for discrete Morse complexes on tetrahedral meshes},
  journal   = {Comput. Graph. Forum},
  volume    = {32},
  number    = {3},
  pages     = {361--370},
  year      = {2013},
  url       = {http://dx.doi.org/10.1111/cgf.12123},
  doi       = {10.1111/cgf.12123}}
```
if you used the triangle based gradient encoding or

```
@article{Comi16,
  author    = {Lidija Comic
               Leila De Floriani and
               Federico Iuricich and
               Paola Magillo},
  title     = {Computing a Forman Gradient From a Watershed decomposition},
  journal   = {Computetr \& Graphics, to appear.}
```
if you used the watershedToForman algorithm
