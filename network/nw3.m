function nw3(n)
%NW3     Investigation of a singular normal equation matrix
%        for a free leveling line

%Kai Borre June 2, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

if nargin == 0
   n = 10;
end

A = diag(2*ones(n,1))+diag(-ones(n-1,1),1)+diag(-ones(n-1,1),-1);
% The upper left and lower right components 
% are modified corresponding to free boundaries
A(1,1) = 1;
A(n,n) = 1;
% We calculate the eigenvectors V and the eigenvalues D
[V,D] = eig(A);
% We set the matrix of eigenvectors Psi
Psi = zeros(n,n);
% and put the correct results into it
Psi(:,1) = ones(n,1)/sqrt(n);
for j = 2:n
   for i = 1:n
      Psi(i,j)=cos((2*i-1)*(j-1)*pi/(2*n))*sqrt(2)/sqrt(n);
   end;
end;

% Finally the eigenvalues are calculated
for j = 1:n
   lambda(j) = 4*(sin(((j-1)*pi)/(2*n)))^2; 
end
% Calculate the pseudoinverse of the original matrix A
Apl = pinv(A);
% The pseudoinverse according to our formula
Aplus = zeros(n,n);

for i = 1:n
   for j = 1:n
      for k = 1:n-1
         Aplus(i,j) = Aplus(i,j)+cos((2*i-1)*k*pi/(2*n))*...
            cos((2*j-1)*k*pi/(2*n))/...
            (2*n*(sin((k*pi)/(2*n)))^2); 
      end;
   end;
end

% We output some variables
if n < 10   
   V
   Psi
   diag(D)
   lambda
   Apl
   Aplus
end

% and some graphics
figure(1);
meshz(Psi)
view(-30,20), pause(3)
figure(gcf)
view(-30,60), pause(3)
figure(gcf)
view(-30,80), pause(3)
print network1 -deps

figure(2);
figure(gcf)
meshz(Psi)
figure(gcf)
meshz(A)
figure(gcf)
meshz(Aplus)
figure(gcf)
view(-160,30), pause(3)
figure(gcf)
%%%%%%%%%%%%%%%%% end nw3.m  %%%%%%%%%%%%%%%%%
