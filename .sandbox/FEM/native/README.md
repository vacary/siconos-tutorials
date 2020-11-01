# An attempt to an light weight implementation of FEM in siconos

## source code ./src

This directory contains soure code of the light wieght implementation of FEM in siconos. Its implement simple isoparametric element (T3 and TH4 for the moment)

### Mesh.* (.hpp .cpp)

A class for storing the geometrical properties of the Mesh (vertices, elements).

### MeshUtils.*

Some function to handle meshes:

- createMeshFromGMSH2 a reader of gmsh v2 file (for other version or mesh data format you can use meshio python module to convert into GMSH v2)
- prepareWriteDisplacementforPython, writeDisplacementforPython output of displacement and mesh to read it in python
  


### FiniteElementLinearTIDS

Finite Element discretization of elastic solids that inherits from Lagrangian Linear Systems with time invariant coefficients
  - $M\dot v + Cv + Kq = F_{ext}(t,z) + p $
 
  
### FiniteElementModel

A class that build elementary matrices from isoparametric FEM.

## Simple Examples

### T3.cpp

### TH4.cpp

