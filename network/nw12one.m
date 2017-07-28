function nw12one(n)
%NW12ONE  Spectral density function for 1-D Laplacian.
%        n must be 100 or more

%Kai Borre June 3, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

if nargin == 0, n = 10000; end

l = zeros(n,1);
for i = 1:n
   l(i) = 4*(sin((i-1)*pi/(2*n)))^2; 
end

t = linspace(0,4,n);
td = linspace(0,4,100);
N = zeros(100,1);
j = 1;
for i = 1:100
   while l(j) <= td(i) & j < n
      N(i) = N(i)+1;
      j = j+1;
   end;
end
tt = linspace(0,4,n-1);
figure;
hold on
plot(t,l,'b')		                       % density distribution curve D
plot(tt,diff(l)' ./diff(t),'g')          % curve for differential of D
plot(td(1,2:100), N(2:100,1)/sqrt(n),'r')% curve for integral of D
hold off
%%%%%%%%%%%%%%%%%%% end nw12one.m  %%%%%%%%%%%%%%%%%%%%%%
