function nw12two(m,n)
%NW12TWO  Spectral density function for 2-D Laplacian

%Kai Borre June 3, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

if nargin == 0, m = 200; n = 200; end

lambda = zeros(m,n);
for j = 1:n,
  for i = 1:m,
   lambda(i,j) = 4*((sin((i-1)*pi/(2*m)))^2+(sin((j-1)*pi/(2*n)))^2);
 end
end
l = reshape(lambda,m*n,1);
y = sort(l);
hist(y,400)
%%%%%%%%%%%%%%%%%% end nw12two.m %%%%%%%%%%%%%%%%%%%%%%
