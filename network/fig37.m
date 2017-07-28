%FIG37   Script for plotting Figure 3.7: Confocal ellipses

%Kai Borre April 17, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

x = 0:.01:4;
G = [];
for u = 0:.5:2
   ell = real(sinh(u)*sqrt(real(1 - x.^2/(cosh(u))^2)));
   G = [G ell'];  
end

plot(x,G,-x,G,x,-G,-x,-G,1,0,'o',-1,0,'o')
axis('off')
text(.95,-.4,'1','FontSize',18)
text(-1.1,-.4,'-1','FontSize',18)
text(2.2,-1,'\itu\rm = 1.5','FontSize',18)
print -deps fig37
%%%%%%%%%%%%%%%%% end fig37.m %%%%%%%%%%%%%