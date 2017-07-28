function [P,Pf] = legen_nm(n,m,t)
%LEGEN_NM   We implement the formula (1-62) in the reference for 
%           computing the associated Legendre function Pnm(t) of
%           degree n and order m.
%   
% Reference:  Weikko A. Heiskanen & Helmut Moritz (1967): Physical 
%             Geodesy, W.H. Freeman

%Written by Kai Borre
%March 25, 2000

r = floor((n-m)/2);
sum = 0;
for k = 0:r
   sum = sum+(-1)^k*(prod(2:2*(n-k))/...
      (prod(2:k)*prod(2:n-k)*prod(2:n-m-2*k)))*t^(n-m-2*k);
end
Pf = 2^(-n)*(1 - t^2)^(m/2)*sum;
k = 2;
if m == 0, k = 1; end
P = Pf*sqrt((prod(2:n-m)*(2*n+1)*k)/(prod(2:n+m)));
%%%%%%%%%%%%%%%%%%%%%%% end legen_nm.m  %%%%%%%%%%%%%%%%