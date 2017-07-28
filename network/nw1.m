%NW1       Example 2 on page 185 in
%               G. Strang: Introduction to Linear Algebra

%Kai Borre March 1, 2000
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

A = [1 0 0;1 1 1;1 2 4];
b = [6;0;0];
x = A\b;
t = linspace(0,2);    % by default t is 100 by 1 vector
y = x(1)+x(2).*t+x(3).*t.^2;
plot(t,y)
%%%%%%%%%%%%%%% end nw1.m  %%%%%%%%%%%%

