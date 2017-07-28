function yy = lw(x)
%LW      Spectral density function yy for the discrete Laplacian 
%        operator in 2-D.

%Kai Borre April 27, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

global T

%yy = asin(sqrt((T/4)*ones(1,size(x))...
%        - (T/4)*x .^2)) ./(sqrt(abs(ones(1,size(x)) ...
%                                - (T/4)*x .^2)+100*eps));
yy = asin((sqrt(1-x.^2)*T/4)./sqrt(1-x.^2*T/4));
%%%%%%%%%%%%%%%%%% end lw.m  %%%%%%%%%%%%%%%%%%%%%%