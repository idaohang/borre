function levrecta(m,n)
%LEVRECTA Plot of standard deviations (sigmas) of heights in an m by n 
%         regular, rectangular, free leveling network

%         Implementation of formulas to be found in 
%             K. Borre & P. Meissl:
%             Strength Analysis of Leveling-Type Networks,
%             An Application of Random Walk Theory, Appendix B
%             Geod\aetisk Institut, Meddelelse No. 50, K\obenhavn 1974

%Kai Borre June 1, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

Sigma = zeros(m,n);
for j = 1:m-1
   alpha(j) = acosh(2-cos(j*pi/m)); 
end

for p = 1:m
   for q = 1:n
      Sq = 0;
      for k = 1:n-1
         Sq = Sq+(cos((q+1.5)*k*pi/n))^2/(2*m*n*(sin(k*pi/(2*n)))^2); 
      end
      for j = 1:m-1
         Sigma(p,q) = Sigma(p,q)+ ...
                  2*(cos((p+1.5)*j*pi/m))^2*cosh((q+1.5)*alpha(j))* ...
                     cosh((n-1.5-q)*alpha(j))/ ...
                      (m*sinh(alpha(j))*sinh((n+1)*alpha(j)));
      end
      Sigma(p,q) = Sigma(p,q)+Sq;
   end
end

for i = 1:m
   for j = 1:n
      sigma(i,j) = sqrt(Sigma(i,j));
   end
end

figure;
subplot(211), cp = contour(sigma); clabel(cp)
title('Standard deviations in 2-D free levelling network')
subplot(212), mesh(sigma)
print levrecta -deps
%%%%%%%%%%% end levrecta.m  %%%%%%%%%%%%%%%%%