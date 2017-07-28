function [L,D] = ldlt(A)
%        Algorithm 4.2.1 in Golub and Van Loan

% Written by Kai Borre
% March 21, 2002

if nargin < 1
    A = [10 20 30;20 45 80;30 80 171];
end

n = size(A,1);
v = zeros(n,1);
for j = 1:n
    for i = 1:j-1
        v(i,1) = A(j,i)*A(i,i);
    end
    v(j,1) = A(j,j)-A(j,1:j-1)*v(1:j-1,1);
    A(j,j) = v(j,1);
    A(j+1:n,j) = (A(j+1:n,j)-A(j+1:n,1:j-1)*v(1:j-1,1))/v(j,1);
end
L = tril(A,-1)+eye(n);
D = diag(A);
%%%%%%%%%%%%%%% end ldlt.m  %%%%%%%%%%%%%%  