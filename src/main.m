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
phi_True = @(t,x,y) (1./(4.*t+1)) ...
             .* exp(-((x-ux.*t-0.5).^2+(y-uy.*t-0.5).^2)./(v.*(4.*t+1)));

% --- CONTOUR PLOT ---
% number of interior levels
L = 5;

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

% recalculate dt now that Nt is known (if dt wasn't divisible)
dt = (tmax - tmin) / (Nt-1);

% spatial arrays
x = linspace(xmin, xmax, Nx);
y = linspace(ymin, ymax, Ny);

% temporal array
t = linspace(tmin, tmax, Nt);

%% SOLVER

% --- ALLOCATE 2D phi ARRAY ---
phi = zeros(Nx,Ny);

% --- APPLY INITIAL CONDITION ---
for j = 1:Ny
    for i = 1:Nx
        phi(i,j) = exp(- ((x(i) - 0.5)^2 + (y(j) - 0.5)^2) / v);
        % or just phi(i,j) = phi_True(0,x(i),y(j));

        % In this particular problem this also handles BC's because they
        % are Dirichlet bounds of the analytical solution 
    end
end

% --- APPLY BOUNDARY CONDITIONS --- 
% (because normally your initial condition doesn't perfectly set up BCs)

% set initial time
n = 1;

% vertical bounds
% left
phi(1,:)    = phi_True(t(n),xmin,y);
% right
phi(Nx,:)   = phi_True(t(n),xmax,y);

% horizontal bounds
% bottom
phi(:,1)    = phi_True(t(n),x,ymin).';
% top
phi(:,Ny)   = phi_True(t(n),x,ymax).';

%  --- EXPLICIT EULER ---
% allocate space for new (n+1) phi
phinew = phi;
% iterate over time
for n = 1:Nt-1 % (Nt-1 because the last time is already computed after)
    % compute the domain
    for j = 2:Ny-1
        for i = 2:Nx-1
            phinew(i,j) = phi(i,j) - dt * ...
                ( ux/(2*dx)*(phi(i+1,j) - phi(i-1,j))   ...
                  + uy/(2*dy)*(phi(i,j+1) - phi(i,j-1)) ...
                  - v/dx^2*(phi(i+1,j) - 2*phi(i,j) + phi(i-1,j)) ...
                  - v/dy^2*(phi(i,j+1) - 2*phi(i,j) + phi(i,j-1)) ...
                );
        end
    end

    % reapply bounds
    % vertical bounds
    % left
    phinew(1,:)    = phi_True(t(n+1),xmin,y);
    % right
    phinew(Nx,:)   = phi_True(t(n+1),xmax,y);
    % horizontal bounds
    % bottom
    phinew(:,1)    = phi_True(t(n+1),x,ymin).';
    % top
    phinew(:,Ny)   = phi_True(t(n+1),x,ymax).';
    
    % reset phi
    phi = phinew;
end

%% ANALYTICAL SOLUTION (at t = tmax)

[X,Y] = ndgrid(x,y); % ndgrid used to conserve i,j grid style
phi_exact = phi_True(tmax,X,Y);


%% RESULTS FIGURE

% --- COMMON CONTOUR SYSTEM ---

% generate common contour levels so numerical and analytical contours
phiMin = min([phi(:); phi_exact(:)]);
phiMax = max([phi(:); phi_exact(:)]);

% ignore bounds because they are outliers
allLevels = linspace(phiMin,phiMax,L+2);
levels = allLevels(2:end-1);

% --- PLOT FIGURE ---
% (transpose phi and phi_exact when plotting due to
%  MATLAB plotting convention)

% create figure
fig = figure;
theme(fig,'light');

% plot numerical contours
[c,h] = contour(x,y,phi.',levels, 'LineStyle','-', 'LineWidth',1.5);
hold on;

% plot analytical contours
[c_true,h_true] = ...
    contour(x,y,phi_exact.',levels, 'LineStyle',':', 'LineWidth',1.5);

% --- FIGURE FORMATTING ---

% color bar
clim([phiMin phiMax]);
cb = colorbar;
cb.Label.String = '\phi';

% contour labels
clabel(c_true,h_true,'FontSize',12, 'LabelSpacing',200);
% only analytical gets labels, otherwise its too cluttered

% axis labels
xlabel('x');
ylabel('y');

% title
title(...
    sprintf('Numerical vs. Analytical Solutions for ϕ (t = %.1f)', tmax));

% legend
legend([h,h_true], {'Numerical','Analytical'}, ...
    'Location','south');

% axis limits
xlim([xmin xmax]);
ylim([ymin ymax]);
axis equal;

% final details
grid on;
box on;
set(gca,'FontSize',11);
hold off;

% increase font size
fontsize(fig,"scale",1.2);

% --- EXPORT FIGURE --- 

exportgraphics(fig,'HW1_Contour_Comparison.png', 'Resolution',600, ...
    'BackgroundColor','white');

%% ERROR ANALYSIS
% vectorized using MATLAB

% total number of spatial nodes
N = Nx * Ny;

% point-wise error
error = phi_exact - phi;

% RMS error
epsilon = sqrt(sum(error(:).^2) / N);

% display RMS error
fprintf('\nRMS Error at t = %.2f: %.8e\n', tmax, epsilon);