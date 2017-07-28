function draw_ell(X,Y,A)
%DRAW_ELL  Plot at (X,Y) of ellipse given by the 2 by 2 matrix A

%Copyright (c) by Kai Borre
%$Revision: 1.1$  $Date:2000/04/21  $

[m,n] = size(A);
if m ~= 2 | n ~= 2 error('Wrong dimension of matrix'); end
[v,d] = eig(A);
if d(1,1) <= 0 | d(2,2) <= 0,
   error('The input matrix is not positive definite'); end;

% Computation for ellipse
[lambda,k] = sort(diag(d));
v = v(:,k);
if any(any(v)) == 1
   alpha = atan2(v(2,2),v(1,2));
else
   alpha=0; 
end
rot = [cos(alpha) sin(alpha);-sin(alpha) cos(alpha)];
t = linspace(0,2*pi,100);
a = sqrt(lambda(2));
b = sqrt(lambda(1));
pl = [a*cos(t);b*sin(t)];
for t = 1:100
   current = rot*pl(:,t);
   curve(1:2,t) = current; 
end

plot(X+curve(2,1:100),Y+curve(1,1:100),'-')
axis('equal')
%%%%%%%% end draw_ell.m %%%%%%%%%%%%%%%%%%%
