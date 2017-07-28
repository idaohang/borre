% Solution to Eksempel 6.9
% October 31, 1999

A = [1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1;
     1 1 0 0;
     0 1 -1 0;
     0 0 1 1];
cov = diag([1. 1.11 .91 1.25 1.0 1.11 1.0]);
b = [20.21;40.07;34.17;35.84;60.40;5.87;69.99];
[x,sigma_x] = lscov(A,b,cov);
[m,n] = size(A);
for i = 1:n, fprintf('\n x(%d): %8.4f m',i,x(i)), end
fprintf('\n')
U = (chol(cov))';
z = U\b;
W = U\A;
mse = (z'*z-x'*W'*z)/(m-n);
fprintf(['\n Standard deviation for unit weight:',...
          ' %6.4f m/sqrt(km)\n'],sqrt(mse))
p = A*x;
for i = 1:m
   fprintf('\n Estimated observations b(%d): %8.4f m',i,p(i))
end
fprintf('\n')
fprintf('\n Test for loop sum: %6.4f m', p(1)-p(5)+p(2))
fprintf('\n Test for loop sum: %6.4f m', -p(2)+p(6)+p(3))
fprintf('\n Test for loop sum: %6.4f m', -p(3)+p(7)-p(4))
fprintf('\n')
db = sqrt(diag(A*inv(A'*inv(cov)*A)*A')*mse);
for i = 1:m, fprintf('\n sigma_obs(%d): %5.1f mm',i,...
                       1000*db(i)), end
fprintf('\n')
%%%%%%%%%%%%%%%%%%%%%  eks69.m  %%%%%%%%%%%%%%%
