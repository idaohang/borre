function [G,Nplus] = nw8(n)  
%NW8       Discrete approximation of continuous 1-D leveling with
%          Neumann boundary conditions

%Kai Borre June 5, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

if nargin == 0, n = 10; end

u = linspace(1/(2*n),1-1/(2*n),n);
v = u;
G = zeros(n,n);
for i = 1:n
   for j = 1:n
      if v(j) > u(i)
         G(j,i) = (u(i)^2+v(j)^2)/2+1/3-v(j);
      else
         G(j,i) = (u(i)^2+v(j)^2)/2+1/3-u(i); 
      end
   end
end

Nplus = zeros(n,n);
for i = 1:n
   for j = 1:n
      if j >= i
         Nplus(i,j) = ...
            ((n-1)*(2*n+5-6*j)/3+(j-1)*(j-2)+i*(i-1))/(2*n);
      else
         Nplus(i,j) = ((n-1)*(2*n+5-6*i)/3 ...
                                +(i-1)*(i-2)+j*(j-1))/(2*n);  
      end
   end
end

if n < 10
   G*n-Nplus    % approximation error
end
plot(G*n,'r')
hold on
pause
plot(Nplus,'go')
hold off
print -deps nw8 
s%%%%%%%%%%%%%%% end nw8.m  %%%%%%%%%%%%