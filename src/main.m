% MATLAB Code for CTP HW 1, 09/01/2026
% Author: Melik Demirel
% Carnegie Mellon University

clear; close all; clc;

%% GIVEN PARAMETERS

% --- GEOMETRY GIVENS ---

% x bounds
xmin = 0; 
xmax = 2;

% y bounds
ymin = 0;
ymax = 2;

% time bounds
tmin = 0;
tmax = 1.0; % (this is also T)

% CFL < 1 (to calculate dt using a convective restriction only)
CFL = 0.01;

% node counts
Nx = 81;
Ny = 81;

% --- PHYSICS GIVENS ---

% convection velocities
ux = 0.8;
uy = 0.8;

% viscosity
v = 0.01; % nu

% --- ACTUAL SOLUTION ---

% function handle
phi_True = @(t,x,y) (1/(4*t+1)) ...
                    * exp(-((x-ux*t-0.5)^2+(y-uy*t-0.5)^2)/(v*(4*t+1)));


%% SPATIAL AND TEMPORAL ARRAYS

% cell counts (intervals calculated from node counts)
xCells = Nx-1;
yCells = Ny-1;

% spatial steps
dx = (xmax - xmin) / xCells; 
dy = (ymax - ymin) / yCells;

% temporal step calculation using convective CFL restriction
dtx = dx * CFL / abs(ux);
dty = dy * CFL / abs(uy);
dt = min(dtx,dty);

% time steps is found by reversing the process for space
Nt = ceil((tmax - tmin) / dt) + 1; % (+1 converts intervals to nodes)
% ceil is used to ensure that in the worst case (when the range isn't 
% divisible well by dt), there are more time-steps than necessary

% spatial arrays
x = linspace(xmin, xmax, Nx);
y = linspace(ymin, ymax, Ny);

% temporal array
t = linspace(tmin, tmax, Nt);


%% SOLVER

% allocate 2D phi array




