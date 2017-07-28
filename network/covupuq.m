%COVUPUQ  Plot of covariance function for 
%         relative direction network cov(u_P,u_Q)

%Kai Borre March 1, 2000
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

l = 1;
sigma2 = 1;

for i = 1:20      % x
   for j = 1:20   % y
      [theta,r] = cart2pol(i,j);
      a1 = sqrt(2)*r/l;
      cov(i,j) = 2*(3*log(r)+besselk(0,a1)- ...
                   cos(2*theta)*(besselk(0,a1)+ ...
                     a1*besselk(1,a1)-(l/r)^2))/(3*sqrt(3)*pi)*sigma2;
   end
end
surfl(cov)
%%%%%%%%%%%%% end covupuq.m  %%%%%%%%%%%%%%%%%%