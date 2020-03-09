# theofun

theofun is a collection of MATLAB functions designed to create theoretical morphospaces from 2D EFA data. Theoretical shapes can then be tested for functional performance, using FEA.

List of functions:

- taxaSetFromFile: generate taxaSet object from raw data.
- taxaSet: main class to hold/manipulate raw EFA and taxon datasets.

- landscape: class that holds x,y and z grid data. Used to generate surface/color map plots
- meshSpace: class that holds a grid of mesh2D objects across a theoretical morphospace, for generation of performance surfaces
- pareto: class that holds a 2D pareto dataset, and performs pareto ranking functions
- Phylogeny: phylogeny class, for plotting (can't be used for any phylogenetic testing)
- shapespace: theoretical morphospace, with grid of theoretical shapes (objects of class theoShapeN)
- timeMorph: time-binned morphospace data, used for disparity analyses and visualisation.

- ele2D (internal class): class to hold element data for 2D meshes
- mesh2D (internal class): class to hold 2D meshes and perform functional analyses
- meshList (internal class): list of mesh2D objects, perform bulk functional analyses
- PhyloNode (internal class): used in construction and functions of Phylogeny objects
- taxonData (internal class): individual taxon data.
- theoShapeN (internal class): individual theoretical shape data.

All other functions are internal, and described within.
