%EX8     Kalman filtering of usual least-squares problem.
%        A is the coefficient matrix, b observations 
%        of equal weight. 
%        With increasing i we filter more observations 

%Kai Borre 25-11-98
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 1998/11/25  $

A = [1 1;1 2;-1 1];
b = [2;1;0];
%Covariance of observations
Cov = diag([1 .5 1]);
%Initial state vector
x = zeros(size(A,2),1);
%Initial weight
P = 1.e6*eye(size(A,2)); 

for i = 1:size(b,1)
   [x,P] = k_update(x,P,A(i,:),b(i),Cov(i,i))
end

fprintf('\nSolution:\n');
   for j = 1:size(A,2)
      fprintf('  x(%2g) = %6.3f\n',j,x(j));
   end
%%%%%%%%%%%%%%%%% end ex8.m %%%%%%%%%%%%%%%%%%%




