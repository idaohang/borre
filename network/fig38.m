%FIG38   Script for plotting Figure 3.8
%        Ellipses tending to fill out an infinite parallel strip

%Kai Borre April 17,1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

x = 0:.02:8;
G = [];
for u0 = 1.0:-.1:.2
   ell = pi*real(sqrt(1-x.^2/(pi*cosh(u0)/(2*sinh(u0)))^2))/2;
   G = [G ell'];  
end

h = axes('PlotBoxAspectRatio',[ 1 .2 1]);
plot(x,G,-x,G,x,-G,-x,-G, ...
          [-8 8],[pi/2 pi/2],'--',[-8 8],[-pi/2 -pi/2],'--')
axis('off','equal')
print -deps fig38
%%%%%%%%%%%%%%%%% end fig38.m %%%%%%%%%%%%%
