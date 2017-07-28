%CLOWN  Demonstration of image compression by means of SVD
%       We illustrate low-rank approximations of a clown

%Copyright (c) by Kai Borre
%$Revision: 1.0$   $Date:1999/10/30  $

load clown.mat;
[U,S,V] = svd(X);
colormap('gray');
for k = [3,10,20,200]
   image(U(:,1:k)*S(1:k,1:k)*V(:,1:k)')
   pause
end
%%%%%%%%%%%%%%%%%%%%%% end clown.m  %%%%%%%%%%