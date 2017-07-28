%FIG51  Plot of 1-dimensional density function for 
%       the discrete Laplacian

%Kai Borre April 27, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

x = 0.01:.005:3.999;
my = 1 ./(pi*sqrt(x.*(4-x)));

plot(x,my)
xlabel('\ity','FontSize',18)
ylabel('\it\mu','FontSize',18,'VerticalAlignment','Bottom')
set(gca,'FontSize',18)
print -deps fig51
%%%%%%%%%%%%%%%% end fig51.m  %%%%%%%%%%%%%%%%%
