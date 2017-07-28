%FIG55  Plots upper and lower bounds for the spectral distribution function
%       N(lambda) for the 2-D density function of the discrete Laplacian

%Kai Borre April 27, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

global T

m = 100;
n = 50;
lambda = 0.01:.025:3.99;
x = 0.0001:.01:.99;
q = size(lambda,2);
N = zeros(1,q);
for i = 1:q
   T = lambda(i);
   Q = quad8('lw',0,1,.05);
   N(1,i) = 2*m*n*sqrt(T)*Q/pi^2;
end

% adding the last half
N = [N 2*N(q)-fliplr(N)];
lambda = [lambda lambda+3.98];

plot(lambda,N,'-')
xlabel('\lambda','Fontsize',18)
ylabel('\itN \rm(\lambda)','FontSize',18,'VerticalAlignment','Bottom')
set(gca,'FontSize',18)
print -deps fig55
%%%%%%%%%%%%%%%%%%%% end fig55.m  %%%%%%%%%%%%%%%%%%%
