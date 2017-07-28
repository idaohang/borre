function nw11(zeta)
%NW11      Green function for a unit disc with Neumann boundary condition.
%          This function has created the cover figure of the present book

%Kai Borre April 26, 1995
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

if nargin == 0, zeta = .5+.6i; end

z = cplxgrid(32);
for j = 1 %:4
   h = figure;
   cplxmap(z,-(log(abs(z-zeta).^2)+log((abs(z-zeta).^2))...
      +(abs(z).^2-1)*(abs(zeta).^2-1)...
      -(abs(z)).^2-(abs(zeta)).^2+3/2)/(4*pi));
   view(j*45,30)
   axis('off')
end

%print -deps greenneu %this command is time consuming
%%%%%%%%%%%%%%%%%%% end nw11.m %%%%%%%%%%%%%%%