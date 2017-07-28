%NW9        Doing all the basic steps in the network theory for a 
%           single triangle. Transformation of observations, computation
%           of covariance matrices and their transformation, etc

%Kai Borre March 1, 2000
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

%Coordinates of nodes
x = [0;1;0];
y = [0;1;1];
dx = diff([x;x(1)])
dy = diff([y;y(1)])
%Observations
b = [1;1;-2.1];
%weights
c(1,1) = 1/norm([x(1)-x(2);y(1)-y(2)]);
c(2,1) = 1/norm([x(2)-x(3);y(2)-y(3)]);
c(3,1) = 1/norm([x(3)-x(1);y(3)-y(1)]);
Sigma_b = diag(1./c)    % (3.2)
V_inv = [dx dy 1./c]
V = inv(V_inv);
f = V*b                 % f(3) is the loop sum!
A = [dx dy]
x_hat = A\(V_inv*f)
Sigma_f = V*Sigma_b*V'
S = inv(Sigma_f)
g = [x_hat;0]-f
E = g'*S*g
r = b-A*x_hat;
E = r'*inv(Sigma_b)*r
%%%%%%%%%%%%%%% end nw9.m  %%%%%%%%%%%%%%%%%%%