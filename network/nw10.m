%NW10     Computing the pseudo-observations for a single triangle
%         when given distance measurements

%Kai Borre March 1, 2000
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

%Coordinates of nodes
x = [0;1;0];
y = [0;1;1];
dx = diff([x;x(1)])
dy = diff([y;y(1)])
%observations
b = [0.350;0.005;0.005];
%pseudo-observations
s(1,1) = norm([x(1)-x(2) y(1)-y(2)]);
s(2,1) = norm([x(2)-x(3) y(2)-y(3)]);
s(3,1) = norm([x(3)-x(1) y(3)-y(1)]);
f = b-log(s)
theta = atan2(dy,dx)
M = [cos(theta(1))^2, sin(theta(1))^2, 2*sin(theta(1))*cos(theta(1));
     cos(theta(2))^2, sin(theta(2))^2, 2*sin(theta(2))*cos(theta(2));
     cos(theta(3))^2, sin(theta(3))^2, 2*sin(theta(3))*cos(theta(3))];
determinant = det(M)
g = M\f
H = [g(1) g(3); g(3) g(2)];
G = eye(2)+H
t = ones(2,1);
scalar = t'*G*t
for i = 1:3
   h(i,1) = [dx(i) dy(i)]*H*[dx(i) dy(i)]'/(s(i))^2;
end
h
%%%%%%%%%%%%%%% end nw10.m  %%%%%%%%%%%%%%%%%%%