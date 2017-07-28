function fig54(c);
%FIG54  Plots upper and lower bounds for the spectral
%       distribution function N(lambda)in case of 
%       1-D density function of the discrete Laplacian 

%Kai Borre April 27, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

n = 100;
x = 0.01:.005:3.999;
N = 2*n*asin(sqrt(x)/2)/pi;
U = N.*(1+c./(n.*sqrt(x)));
L = N.*(1-c./(n.*sqrt(x)));
plot(x,N,'-',x,U,'--',x,L,'--')
xlabel('\lambda','Fontsize',18)
ylabel('\itN \rm(\lambda)','Fontsize',18,'VerticalAlignment','Bottom')
set(gca,'Fontsize',18)
print -deps fig54
%%%%%%%%%%%%%% end fig54.m  %%%%%%%%%%%%%%%%%%%%%%%%