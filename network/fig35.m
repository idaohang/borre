%FIG35   Plot of Figure 3.5
%        Green's function for a unit circle, theta = pi/2

%Kai Borre April 16, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

r = 0:.001:1;
G = [];
for r0 = eps:.1:1+eps
   g = -(log(r0.^2+r.^2)+log(r.^2.*r0.^2+1)-r.^2-r0.^2+3/2)/(4*pi);
   G = [G g'];  
end

plot(r,G)
axis([0 1 -0.1 0.6])
ylabel('\itG','FontSize',18)
xlabel('\itr','FontSize',18)
text(.1,.3,'\itr\rm_0 = 0.0','rotation',-70) 
text(.02,.16,'\itr\rm_0 = 0.2','rotation',-20)
text(.02,-.07,'\itr\rm_0 = 1.0') 
set(gca,'FontSize',18)

print -deps fig35
%%%%%%%%%%%%%%%%% end fig35.m %%%%%%%%%%%%%