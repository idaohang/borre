%NW2    Problem 4.3.1 in
%            G. Strang: Introduction to Linear Algebra

%Kai Borre March 1, 2000
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

A = [1 0;1 1;1 3;1 4];
b = [0;8;8;20];
x = A\b;
p = [1;5;13;17];
y = A\p;
t = [0;1;3;4];

figure
hold on
plot(t,b,'rx')
plot(t,p,'o')
plot([0 4],[x(1) x(1)+x(2)*4],'r-.')
plot([0 4],[y(1) y(1)+y(2)*4],':')
hold off
%%%%%%%%%%%%%%% end nw2.m  %%%%%%%%%%%%

