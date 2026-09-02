% MATLAB Code for CTP HW 1, 09/01/2026
% Author: Melik Demirel
% Carnegie Mellon University

clear; close all; clc;

%% GIVEN PARAMETERS

% GEOMETRY GIVENS

% x bounds
xmin = 0; 
xmax = 2;

% y bounds
ymin = 0;
ymax = 2;

% time bounds
tmin = 0;
tmax = 1.0; % (this is also T)

% time step (chosen, not given)
dt = 0.001;

% node counts
Nx = 81;
Ny = 81;

% PHYSICS GIVENS

% convection velocities
ux = 0.8;
uy = 0.8;

% viscosity
v = 0.01; % nu

%% SIMPLE CALCULATIONS

% cell counts (intervals calculated from node counts)
xCells = Nx-1;
yCells = Ny-1;

% geometric steps
dx = (xmax - xmin) / xCells; 
dy = (ymax - ymin) / yCells;

% geometry arrays
x = linspace(xmin, xmax, Nx);
y = linspace(ymin, ymax, Ny);

% grid
[X,Y] = meshgrid(x,y);

