% GEONET
%   This script solves the system of PDEs for geodetic networks (4.23)
%   with homogeneous boundary conditions (4.29). 
%   The components of the right-hand side of (4.23) are all taken to be
%   delta - 1/(area of region).
%   The right-hand side of (4.29) is 0.
%   The solution components (u, v) are plotted as vector field arrows.
%   The components alpha and beta are plotted as graphs, also colored
%   according to function value (figures 2 and 3, respectively).
%
%   The script uses functions from the finite element software
%   FEMLAB that is distributed by COMSOL.

%Daniel Bertilsson October, 2000
%Copyright (c) by Comsol
%$Revision: 1.0 $  $Date: 2000/12/16 $

%   The following parameters define the problem:
c1 = 1e-6;
c2 = 10;
c3 = 1e-6;
c4 = 1e-6;
l = 0.1;           
% The solution will have oscillations of wavelength comparable to l.
% Thus, if you take a smaller l, you need to refine the mesh in order
% to resolve the oscillations. This is done by changing the value hmax
% below to a smaller value.                  
geo = 'rectangle';   % Geometry type. Can be 'ellipse' or 'rectangle'.
cx = 0.2; cy = 0;    % Coordinates of the center of the ellipse/rectangle. 
lx = 1;   ly = 1;    % Width and height of the ellipse/rectangle. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define variables.
clear fem;
fem.dim = { 'u' 'v' 'alpha' 'beta' };
fem.variables = { 'A' 3*(c1+c2)+c3+c4 ... 
                  'B' c1+c2+3*(c3+c4) ...
                  'a' 4*c2 'b' 4*c4 'K' l^2/3 };
% Define geometry.   
switch geo
 case 'rectangle'
  fem.geom = rect2(cx-lx/2, cx+lx/2, cy-ly/2, cy+ly/2);
  area = lx*ly;
 case 'ellipse'
  fem.geom = ellip2(cx, cy, lx/2, ly/2);
  area = pi*lx*ly/4;  
end
fem.geom = geomcsg( {fem.geom}, {}, {point2(0,0)} );
% Create mesh.
hmax = 0.02;          % Maximum side for the triangles in the mesh.
hmax_orig = 0.0001;   % Local hmax around the origin.
lr = get(fem.geom, 'lr');
nr_orig = find(isfinite(lr{1}));   % Point number in the mesh for the origin.
fem.mesh = meshinit(fem, 'hmax', {hmax, [nr_orig; hmax_orig]}); 
% Natural boundary conditions are default, so we don't need to specify them.
% Define coefficients in the PDE.
fem.equ.c = { { {'A' 'B'}...
                {0 '(A-B)/2'; '(A-B)/2' 0}  {'B' 'A'}...
                0                           0         'K*a'...
                0                           0         0      'K*b' } } ;     
fem.equ.al = { { 0 ...
                 0        0 ...
                 {'a' 0}  {0 'a'}   0 ... 
                 {0 '-b'} {'b' 0}   0   0 } };            
fem.equ.a = { { 0 0 '-2*a' '-2*b' } };
f = -1/area;        
fem.equ.f = { { f f f f } };      
% The delta functions on the right-hand side are specified directly in the
% right-hand side vector L occurring in the FEM discretization of the problem.
np = size(fem.mesh.p, 2);          % Number of points in mesh.
fem.mat.L = zeros(4*np, 1);        
fem.mat.L(nr_orig) = 1;            % The delta function gives a contribution at 
fem.mat.L(nr_orig + np) = 1;       % the origin, for all four components.
fem.mat.L(nr_orig + 2*np) = 1;
fem.mat.L(nr_orig + 3*np) = 1;
% Solve problem
fem.sol = femlin(fem);
% Plot the four components of the solution.
figure
postplot(fem, 'arrowdata', {'u','v'}, 'arrowxspacing', 30, 'arrowyspacing', 30);
figure
postplot(fem, 'tridata', 'alpha', 'triz', 'alpha', 'trimap', 'jet(4096)');
figure
postplot(fem, 'tridata', 'beta', 'triz', 'beta', 'trimap', 'jet(4096)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%  end geonet.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

