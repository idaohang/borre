%FIG36   Script for plotting Figure 3.6
%        Green's function for a unit circle, theta = 0

%Kai Borre April 17, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

r = 0:.001:1;
G = [];
for r0 = 0:.1:1
   g = -(log(r0.^2+r.^2-2.*r.*r0+eps)+log(r.^2.*r0.^2+1-2.*r.*r0+eps)...
           -r.^2-r0.^2+3/2)/(4*pi);
   G = [G g'];  
end

plot(r,G)
axis([0 1 -0.1 0.6])
ylabel('\itG','Fontsize',18)
xlabel('\itr','Fontsize',18)
%right curve branches
text(.9,.53,'\itr\rm_0 = 0.8','Rotation',-30)
text(.9,.37,'\itr\rm_0 = 0.7','rotation',-15)
text(.9,-.06,'\itr\rm_0 = 0.0','Rotation',-2)
%left curve branches
text(.08,.4,'\itr\rm_0 = 0.1','rotation',78)
text(.17,.36,'\itr\rm_0 = 0.2','rotation',75)
text(.42,.12,'\itr\rm_0 = 1.0','rotation',34)
set(gca,'FontSize',18)
print -deps fig36
s%%%%%%%%%%%%%%%%% end fig36.m %%%%%%%%%%%%%