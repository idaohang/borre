function rec_baye(A,b,Sigma)
%REC_BAYE Recursive Least Squares solved via a Bayes Filter.
%         A is the coefficient matrix, b the observations and 
%         Sigma a vector containing the diagonal entries of
%         the covariance matrix for the problem.

%Written by Kai Borre
%August 21, 1999

if nargin == 0
   A = [1 0;1 1;1 3;1 4];
   b = [0;8;8;20];
   Sigma = diag([1,1,1,1]);
end

% Initialization
x = ones(size(A,2),1);
P = 1.e6*eye(size(A,2));

for i = 1:size(b,1)
   [x,P]= b_row(x,P,A(i,:),b(i),Sigma(i,i));
   fprintf('\nSolution:\n');
   for j = 1:size(A,2)
      fprintf('  x(%2g) = %6.3f\n',j,x(j));
   end
end

dof = size(b,1)-size(A,2);
if dof ~= 0
   P = (norm(b-A*x))^2*P/dof;
else
   P = (norm(b-A*x))^2*pinv(A'*Sigma*A);  
end      
fprintf('\nFinal Covariance matrix:\n');
for j = 1:size(A,2)
   for k = 1:size(A,2)
      fprintf('%12.7f',P(j,k));
   end
   fprintf('\n');
end
fprintf('\nTrace of Covariance matrix: %12.7f\n',trace(P));
%%%%%%%%%%%%%%%%% end rec_baye.m %%%%%%%%%%%%%%%%%%%




