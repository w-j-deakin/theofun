function [output] = TemplateFunctionalCalculation(mesh,var)

%% THIS IS A TEMPLATE FUNCTION, FOR USE IN meshSpace.customFunction. %%
% copy this file, edit, rename and save, then call it in desired function.

% INPUT 1: mesh. mesh2D object (see mesh2D.m). This contaons all shape data
% you could need for each theoretical shape.

% INPUT 2: var. user defined. any custome inputs the user would like for
% their function should be input here.

% OUTPUT: user defined. can be single value, array of many values, object.
% Will be stored in a cell array by meshSpace.customFunction


%% EXAMPLE: calculate mass of shape, based on the mesh and an input density
% Note - you should delete this section and replace with your own custom
% code

% calculate area
A = polyarea(mesh.x(mesh.outline), mesh.y(mesh.outline));

% mass = area * density
output = A * var;


end

