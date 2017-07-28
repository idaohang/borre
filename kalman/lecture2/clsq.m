function [x_q,x_s,r_s,x,r,diff] = clsq(A,b,B,d)
%CLSQ     Constrained least squares.
%         Given an overdetermined system Ax = b, rank(A) = n, 
%         with constraints Bx = d, rank(B) = p. 
%         Solution via the generalized Singular Value 
%         Decomposition:  A = UCX'
%                         B = VSX'
%                         C'*C+S'*S = I
% 
%         See
%         Gene H. Golub & Charles F. van Loan:
%                Matrix Computations, Second edition,
%                Section 12.1.1, and
%         \AA{}ke Bj\"orck (1996): Numerical Methods for 
%                Least Squares Problems. SIAM 


%Written by Kai Borre
%November 16, 1998

if nargin == 0
   A = [1 0; 1 1; 1 3];
   b = [0; 8; 8];
   B = [1 4];
   d = [20];
end
m = size(A,1);
n = size(A,2);
p = size(B,1);

%Solution by means of a QR decomposition
[Q,R] = qr(B');
y = R(1:p,1:p)\d;
AQ = A*Q;
z = AQ(:,p+1:n)\(b-AQ(:,1:p)*y);
x_q = Q(:,1:p)*y+Q(:,p+1:n)*z;

%Solution by means of generalized
%singular value decomposition
[U,V,X,C,S] = gsvd(A,B);
q = size(find(diag(C)==0),1);
W = inv(X');
x_s = zeros(n,1); % m
for i = 1:p
   x_s = x_s+V(:,i)'*d/S(i,i)*W(:,i);
end
for i = p+1:n
   x_s = x_s+U(:,i)'*b/C(i,i)*W(:,i);
end
r_s = b-A*x_s;

%For reasons of comparison we quote the 
%solution without constaints
x = zeros(n,1);
for i = q+1:p
   x = x+U(:,i)'*b/C(i,i)*W(:,i);
end
for i = p+1:n
   x = x+U(:,i)'*b*W(:,i);
end
r = b-A*x;

%Sum of constrained - unconstrained residuals squared
diff = 0;
for i = q+1:p
   diff = diff+(U(:,i)'*b-C(i,i)/S(i,i)*V(:,i)'*d)^2;
end
%%%%%%%%%%%%%%%% end clsq.m  %%%%%%%%%%%%%%%%%%%%



