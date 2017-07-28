function x = wlsq(A,b,B,d)
%WLSQ   Constraining a least squares problem by big weight 
%       of observation(s).
%       Given the system Ax=b, rank(A)=n, 
%       and additionally observations with big weight
%       Bx=d, rank(B)=p. 
%       See
%          Charles L. Lawson & Richard J. Hanson:
%          Solving Least Squares Problems, Prentice-Hall, 
%          Inc., 1974, Chapter 22

%Written by Kai Borre
%November 16, 1998

if nargin == 0
   A = [1 0;1 1;1 3];
   b = [0;8;8];
   B = [1 4];
   d = [20];
end

big = 1.e10;
x = (A'*A+big*B'*B)\(A'*b+big*B'*d);
%%%%%%%%%%%%%%%%% end wlsq.m  %%%%%%%%%%%%%%%%%%%%



