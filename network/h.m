function res = H(n)
%H       Generates an H matrix that describes a 2-D leveling network

%Kai Borre June 1, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

H = -eye(n)+diag(ones(n-1,1),1);
col = [zeros(n-1,1); 1];
res = [H col];
%%%%%%%%%% end h.m  %%%%%%%%
