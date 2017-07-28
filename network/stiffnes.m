function S = stiffnes(x1, y1, x2, y2, x3, y3)
%STIFFNES  Computation of stiffness matrix for the triangle
%          with vertices (x1,y1),(x2,y2), and (x3,y3) 

%Kai Borre August 12, 2000
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

% Default weights of the six observations
c = [1; 1; 1; 1; 1; 1];

if nargin == 0
x = [0; 0; 1];
y = [0; 1; 0];
else
x = [x1; x2; x3];
y = [y1; y2; y3];
end
x0 = mean(x);
y0 = mean(y);

for edge = 1:6
if edge == 1, i = 1; j = 2; end; 
if edge == 2, i = 2; j = 3; end;
if edge == 3, i = 3; j = 1; end;
if edge == 4, i = 2; j = 1; end;
if edge == 5, i = 3; j = 2; end;
if edge == 6, i = 1; j = 3; end;
theta(edge,1) = atan2(y(i)-y(j),x(i)-x(j));
theta = theta+pi;
S2T(edge,1) = sin(2*theta(edge));
C2T(edge,1) = cos(2*theta(edge));
end

x = [x; x(2); x(3); x(1)];
y = [y; y(2); y(3); y(1)];
deltax = x-x0;
deltay = y-y0;

S = zeros(6,6);
S(1,1) = c'*c;
S(2,1) = c'*C2T;
S(3,1) = c'*S2T;
S(2,2) = c'*(C2T.^2);
S(3,2) = c'*(C2T.*S2T);
S(4,2) = c'*(deltax.*C2T);
S(5,2) = c'*(deltay.*C2T);
S(3,3) = c'*(S2T.^2);
S(4,3) = c'*(deltax.*S2T);
S(5,3) = c'*(deltay.*S2T);
S(4,4) = c'*(deltax.^2);
S(5,4) = c'*(deltax.*deltay);
S(5,5) = c'*(deltay.^2);
S(6,6) = c'*(1./c);
S = S+tril(S,-1)'; % making S symmetric
%%%%%%%%%%%%%%%%% end stiffnes.m  %%%%%%%%%%%%%

