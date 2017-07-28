function fig16_18(n);
%FIG16_18  Discrete approximation of continuous 1-D leveling with
%          combinations of Dirichlet and Neumann boundary conditions.
%          Number of nodes is n

%Kai Borre March 2, 2000
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

if nargin == 0, n = 10; end

u = 1:n; 
v = u;
G = zeros(n,n);
%fixed-fixed
for i = 1:n
   for j = 1:n
      if v(j) <= u(i)
         G(i,j) = v(j)-u(i)*v(j)/(n+1);
      else
         G(i,j) = u(i)-u(i)*v(j)/(n+1); 
      end
   end
end
G = [zeros(n,1) G zeros(n,1)];
N = zeros(n,n);
for i = 1:n
   for j = 1:n
      if j >= i
         N(i,j) = i*(n+1-j)/(n+1);
      else
         N(i,j) = j*(n+1-i)/(n+1);  
      end
   end
end
figure(1);
hp1 = plot(0:n+1,G,'r');
axis equal
axis([0 n+1 0 2*sqrt(n)])
if nargin == 0 
   axis([0 11 0 3])
end
hold on
plot(N,'go')
hold off
set(gca,'FontSize',18)
set(hp1,'LineWidth',1)
print -deps fig17

% free-fixed
u = 1:n+1;
v = u;
G = zeros(n+1,n+1);
for i = 1:n+1
   for j = 1:n+1
      if v(j) <= u(i)
         G(i,j) = v(j);
      else
         G(i,j) = u(i); 
      end
   end
end

for i = 1:n
   for j = 1:n
      if i < j
         N(i,j) = i;  
      else
         N(i,j) = j;
      end
   end
end
figure(2);
hp2 = plot(G,'r');
hold on
plot([0 1],[0 1],'r')
plot([n n+1],[n n+1],'r')
hold off
axis equal
axis([0 n+1 0 n+1])
hold on
plot(N,'go')
hold off
set(gca,'FontSize',18)
set(hp2,'LineWidth',1)
print -deps fig16

% free-free
G = zeros(n+1,n);
u = linspace(1/(2*n),1-1/(2*n),n);
u = [u 1+1/(2*n)];
v = u;
for i = 1:n+1
   for j = 1:n+1
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
figure(3);
hp3 = plot(G(1:n+1,1:n)*n,'r');
hold on
hp4 = plot([0 1], [flipud(G(1:n,11)*n) G(1:n,1)*n],'r');
%hold off
axis equal
axis([0 n+1 -2*sqrt(n) 2*sqrt(n)])
if nargin == 0
   axis([0 11 -2 3])
end
%hold on
plot(Nplus,'go')
hold off
set(gca,'FontSize',18)
set(hp3,'LineWidth',1)
set(hp4,'LineWidth',1)
print -deps fig18
%%%%%%%%%%%%%%% end fig16_18.m  %%%%%%%%%%%%