# theofun

theofun is a collection of MATLAB functions designed to create theoretical morphospaces from 2D EFA data. Theoretical shapes can then be tested functionally, using FEA.

List of functions:

- ele2D (internal class): class to hold element data for 2D meshes
- landscape: class that holds x,y and z grid data. Used to generate surface/color map plots
- mesh2D (internal class): class to hold 2D meshes and perform functional analyses
- meshList (internal class): List of mesh2D objects, perform bulk functional analyses
- meshSpace: class that holds a grid of mesh2D objects across a theoretical morphospace, for generation of performance surfaces
