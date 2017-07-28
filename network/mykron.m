function [N,O] = mykron(m,n)
%MYKRON     Computation of normals for 2-D leveling network by means of
%           Kronecker products

%Kai Borre June 1, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

if nargin == 0
   m = 3;
   n = 2;
end

if m == 1
   t = m;
   m = n;
   n = t;
end

if n == 1
   if m == 1
      disp('No least-squares problem');
      break
   end
   N = diag(2*ones(m,1))+diag(-ones(m-1,1),1)+ ...
      diag(-ones(m-1,1),-1);
   N(1,1) = 1;
   N(m,m) = 1;
   O = N;
else
% N corresponds to arranging the unknowns  
% row-wise from top to bottom
N = kron(eye(m),h(n-1)'*h(n-1))+kron(h(m-1)'*h(m-1),eye(n));
% O corresponds to arranging the unknowns
% column-wise from left to right
O = kron(h(n-1)'*h(n-1),eye(m))+kron(eye(n),h(m-1)'*h(m-1));
end

pinv(N)
pinv(O)

figure
meshz(pinv(N))
%contourf(pinv(N)), hold on, shading flat
view(30,60)

figure
meshz(pinv(O))
view(30,40)
%%%%%%%%%%%% end mykron.m  %%%%%%%%%%%%
