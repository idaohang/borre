%FIG44   Script for plotting Figure 4.4
%        Bessel functions J_0 and J_1 of first kind

%Kai Borre Ocotber 1, 2000
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

x = 0:.01:10;
plot(x,besselj(0,x))
hold on
plot(x,besselj(1,x))
hold off
axis([0 10 -1 1])
text(1.8,0,'\itJ\rm_0','FontSize',18)
text(3,.4,'\itJ\rm_1','FontSize',18)
set(gca,'FontSize',18)
print -deps fig44
%%%%%%%%%%%%%%%%% end fig44.m %%%%%%%%%%%%%