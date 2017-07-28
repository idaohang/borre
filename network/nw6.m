function nw6(n)
%NW6      Propagation of systematic errors in
%         1-D leveling line with fixed terminals

%Kai Borre February 6, 1995
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

if nargin == 0, n = 10; end
N = diag(2*ones(n,1))-diag(ones(n-1,1),1)-diag(ones(n-1,1),-1);
b = zeros(n,1);
b(ceil(n/2),1) = b(ceil(n/2),1)+1; % introduce error of 1 unit
x0 = inv(N)*b;
%introduce error of .1 unit in a third of observations
b(ceil(n/3):ceil((2*n)/3),1) = b(ceil(n/3):ceil((2*n)/3),1)+.1;
x = inv(N)*b;
plot(x-x0)
%%%%%%%%%%%%%%%%%%%%%%% end nw6.m  %%%%%%%%%%%%%%
