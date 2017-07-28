function ex23(n)
%EX23    Computation of eigensolution of R_N(2,-2,-2)

%Kai Borre June 8, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

global N
N = n;
if nargin == 0, N = 4; end
if N < 3, error('Argument too small'), break, end;

% We input an N x N matrix
A = diag(2*ones(N,1))+diag(-ones(N-1,1),1)+diag(-ones(N-1,1),-1);
% We modify the upper left and lower right components to build
% R_N(2,-2,-2)
A(N,N) = A(N,N)-2;
A(1,1) = A(1,1)-2;
[V,D] = eig(A);
[phi,ind] = sort(diag(D));
x = zeros(N,1);
x(1) = fzero('ex23a',-.5);
x(2) = fzero('ex23b',-.5);
p = ceil(N/2)-1;
for i = 1:p 
   %The following function call needs some refinements
   x(i+2) = fzero('ex23c',i*(2*pi-.9)/(2*N));
end
for i = 1:p
   %The following function call needs some refinements
   x(i+2+p) = fzero('ex23d',i*(2*pi+0.001)/(2*N)); %+0.05
end
lambda = zeros(N,1);
for i = 1:2
   lambda(i) = -4*(sinh(x(i)))^2;
end
for i = 3:N
   lambda(i) = 4*(sin(x(i)))^2; 
end

fprintf('\n');
for i =1:N
   fprintf(' %4.2f',lambda(i));
end
fprintf('\n');
for i = 1:N
   fprintf('\n');
   for j = 1:N
      fprintf('%4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n', ...
                 V(j,ind));
   end
end
fprintf('\n')
clear global
%%%%%%%%%%%%%%%%%%%%%% end ex23.m  %%%%%%%%%%%%%%%%%%%