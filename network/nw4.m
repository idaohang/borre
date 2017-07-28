function nw4(n)
%NW4     Computation and plotting of eigenvectors of
%        the oscillating tri-diagonal matrix (-1,2,-1)

%Kai Borre May 31, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

if nargin == 0
   n = 32;
end

A = diag(2*ones(n,1))-diag(ones(n-1,1),1)-diag(ones(n-1,1),-1);
[v,d] = eig(A);
[y,ind] = sort(diag(d));

meshz(v(:,ind))
view(110,70)
figure;
surf(v(:,ind))
view(80,80)
figure;
waterfall(v(:,ind))
view(20,60)
%contourf(v(:,ind))
s%%%%%%%%%%%%%%%%%%% end nw4.m  %%%%%%%%%%%%