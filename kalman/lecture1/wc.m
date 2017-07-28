% Demonstration of the impact of changing the weights
% in Example 11.11. See also the text of Example 17.9

%Kai Borre 27-6-1997
%Copyright (c) 1997, 1998 by Kai Borre
%$Revision: 1.1 $  $Date: 1998/01/12 $

big = 1.e12;
A = [1 1; 1 2; -1 1];
b = [2; 1; 0];
Cov = diag([1 1/2 1]); 

% Uncomment one of the following lines to see the effect of 
% setting delta C_2 = -1 or delta C_2 = \infty
%Cov(3,3) = big;
%Cov(3,3) = 0;   

x = zeros(2,1);
P = big*eye(2);
for i = 1:3
   [x,P] = k_update(x,P,A(i,:),b(i),Cov(i,i))
   pause
end   

% the same in Bayes' version
x = zeros(2,1);
P = big*eye(2);
for i = 1:3
   [x,P] = b_row(x,P,A(i,:),b(i),Cov(i,i))
   pause
end      
%%%%%%%%%%%%%%% end wc.m  %%%%%%%%%%%%%%%%%%%%%%   