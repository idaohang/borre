function nw7(n)
%NW7      Regular tridiagonal matrices for which we
%         modify the (1,1) and (n,n) components

%Kai Borre June 2, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

echo on
if nargin == 0, n = 4; end

% We input an n by n matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% N_1 = R_n(2,0,0)  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = diag(2*ones(n,1))+diag(-ones(n-1,1),1)+diag(-ones(n-1,1),-1);
% We calculate eigenvectors V and eigenvalues of matrix A
[V,D] = eig(A);
% We define a new matrix for the eigenvectors Psi1
Psi1 = zeros(n,n);
% and put the correct results
for j = 1:n
   for i = 1:n 
      Psi1(i,j) = sqrt(2)*sin(i*j*pi/(n+1))/sqrt(n+1);
   end;
end;
% Finally the eigenvalues are calculated
for j = 1:n 
   lambda(j)=4*(sin((j*pi)/(2*(n+1))))^2; 
end
if n < 10
   V
   Psi1
   L = diag(D);
   L'
   lambda
   pause
end
% Modify the upper left and lower right components to obtain
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% N_2 = R_n(2,-1,0) %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(n,n) = A(n,n)-1;
% We calculate the eigenvectors V and the eigenvalues D
[V,D] = eig(A);
% We define a new matrix for the eigenvectors Psi2
Psi2 = zeros(n,n);
% and put the correct results
% See W.-D. Schuh (1984) Analyse und Konvergenzbeschleunigung der
%	   Methode der konjugierten Gradienten bei geod\"atischen
%	   Netzen. Mitteilungen der geod\"atischen Institute der
%	   Technischen Universit\"at Graz, Folge 49, page 48
for j = 1:n
   for i = 1:n 
      Psi2(i,j)=2/sqrt(2*n+1)*sin((2*j-1)*i*pi/(2*n+1)); 
   end;
end;
% The eigenvalues are calculated, also cf. Schuh
for j = 1:n 
   lambda(j)=4*(sin(((2*j-1)*pi)/(2*(2*n+1))))^2; 
end
if n < 10
   V
   Psi2
   L = diag(D);
   L'
   lambda
   pause
end
% We also modify the upper left components to obtain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% N_3 = R_n(2,-1,-1) %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(1,1) = A(1,1)-1;
% We calculate eigenvectors V and eigenvalues 
[V,D] = eig(A);
% We define a new matrix for the eigenvectors Psi3
Psi3 = zeros(n,n);
% and put the correct results
Psi3(:,1) = ones(n,1)/sqrt(n);
for i = 1:n
   for j = 2:n 
      Psi3(i,j) = sqrt(2/n)*cos((j-1)*(2*i-1)*pi/(2*n));
   end;
end;

% Finally the eigenvalues are calculated
for j = 1:n 
   lambda(j)=4*(sin(((j-1)*pi)/(2*n)))^2;
end
if n < 10
   V
   Psi3
   L = diag(D);
   L'
   lambda
   pause
end
% Calculate the pseudoinverse of the original matrix A
Apl = pinv(A);
% The pseudoinverse according to our special formula
Aplus = zeros(n,n);
% We assign the upper triangular components
for i = 1:n
   for j = i:n
      Aplus(i,j)=((n-1)*(2*n+5-6*j)/3+(j-1)*(j-2)...
         +i*(i-1))/(2*n); 
   end;
end
% Finally we assign the lower triangular components
for i = 1:n
   for j = 1:i
      Aplus(i,j)=Aplus(j,i); 
   end;
end
% We output some interesting numerical results
if n < 10
   Apl
   Aplus
end

h1 = figure;
figure(h1)
plot(Aplus)
title('The Pseudoinverse of {\itA}')

h2 = figure;
figure(h2)
subplot(3,1,1)
meshz(Psi1)
title('Inverses of {\itN_1, N_2,} and {\itN_3}')

subplot(3,1,2)
meshz(Psi2)
subplot(3,1,3)
meshz(Psi3)
view(-30,70)

print -deps nw7
%%%%%%%%%%%%%%% end nw7.m %%%%%%%%%%%%%%