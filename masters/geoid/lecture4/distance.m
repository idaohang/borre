function D = distance(k,l,ml,nl,Deltax,Deltay)
%DISTANCE Computes the kernel for Stoke's integral.
%         The kernel is inversely proportional to 
%         the distance between (k,l) and any grid
%         point. For the grid point (k,l), we set 
%         distance equal to zero.

% Written by Kai Borre
% April, 22, 2000

D = zeros(ml,nl);
for i = 1:ml
   for j = 1:nl
      d = sqrt((i-k)^2*Deltax^2+(j-l)^2*Deltay^2);
      if d == 0, D(i,j) = 0; else D(i,j) = 1/d; end
   end
end
%%%%%%%%%%%%%%%%%%%%% end distance.m  %%%%%%%%%%%%%%%%%%
