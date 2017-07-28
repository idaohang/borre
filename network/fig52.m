%FIG52  Plot of 2-D density function of the discrete Laplacian

%Kai Borre April 27, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

z = 0.01:.005:7.999;
m = sqrt(z .*(8-z)/16);
my = ellipke(m) ./(2*pi^2);
plot(z,my)
xlabel('\itz','FontSize',18)
ylabel('\it\mu','FontSize',18,'VerticalAlignment','Bottom')
set(gca,'FontSize',18)
print -deps fig52
%%%%%%%%%%%%%%%%%% end fig52.m  %%%%%%%%%%%%%%%%%