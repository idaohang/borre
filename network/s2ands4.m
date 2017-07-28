function S = s2ands4(x1, y1, x2, y2, x3, y3, weight2, weight4)
%S2ANDS4  Computation of stiffness matrix for combined
%         type 2 and type 4 observations for the triangle
%         with vertices (x1,y1),(x2,y2), and (x3,y3) 

%Kai Borre August 13, 2000
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

if nargin == 0
x = [0; 0; 1];
y = [0; 1; 0];
weight2 = 1;
weight4 = 1;
else
x = [x1; x2; x3];
y = [y1; y2; y3];
end
x0 = mean(x);
y0 = mean(y);

% Default weights of the observations
c = weight2*ones(6,1);
d = weight4*ones(6,1);

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

S = zeros(8,8);
% We define the lower triangular part of S ~= 0
S(1,1) = c'*ones(6,1);
S(2,1) = c'*C2T;
S(3,1) = c'*S2T;
S(5,1) = -2*c'*deltax;
S(6,1) = -2*c'*deltay;
S(2,2) = c'*(C2T.^2)+d'*(S2T.^2);
S(3,2) = c'*(C2T.*S2T)-d'*(C2T.*S2T);
S(4,2) = d'*S2T;
S(5,2) = -2*c'*(deltax.*C2T);
S(6,2) = -2*c'*(deltay.*C2T);
S(7,2) = -2*d'*(deltax.*S2T);
S(8,2) = -2*d'*(deltax.*S2T);
S(3,3) = c'*(S2T.^2)+d'*(C2T.^2);
S(4,3) = d'*C2T;
S(5,3) = -2*c'*(deltax.*S2T);
S(6,3) = -2*c'*(deltay.*S2T);
S(7,3) = -2*d'*(deltax.*C2T);
S(8,3) = -2*d'*(deltay.*C2T);
S(4,4) = d'*ones(6,1);
S(7,4) = -2*d'*deltax;
S(8,4) = -2*d'*deltay;
S(5,5) = 4*c'*(deltax.^2);
S(6,5) = 4*c'*(deltax.*deltay);
S(6,6) = 4*c'*(deltay.^2);
S(7,7) = 4*d'*(deltax.^2);
S(8,7) = 4*d'*(deltax.*deltay);
S(8,8) = 4*d'*(deltay.^2);
S = S+tril(S,-1)'; % making S symmetric
%%%%%%%%%%%%%%%%% end s2ands4.m  %%%%%%%%%%%%%

