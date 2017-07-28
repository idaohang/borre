%FIG53  Plot of 1-D density function of the discrete bi-harmonic operator 

%Kai Borre April 27, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

y = 0.01:.005:3.999;
my = 1 ./(2*pi^4 .*sqrt(y .^3 .*(y+4-4*sqrt(y))));
plot(y,my)
xlabel('\ity','FontSize',18)
ylabel('\it\mu','Fontsize',18,'VerticalAlignment','Bottom')
set(gca,'Fontsize',18)
print -deps fig53
%%%%%%%%%%%%%%%%%%% end fig53.m  %%%%%%%%%%%%%
