% Demonstration of contour and mesh commands
% from The Matlab 5 Handbook

% Written by Kai Borre
% October 4, 2000

x = 0:.2:3*pi;
y = 0:.25:5*pi;
[XX,YY] =meshgrid(x,y);
Z1 = sin(XX).*sin(YY);

x = -3:.25:3;
y = x;
[XX,YY] = meshgrid(x,y);
Z2 = XX -.5*XX.^3+.2*YY.^2+1;

x = -8:.5:8;
y = x;
[XX,YY] = meshgrid(x,y);
r = sqrt(XX.^2+YY.^2)+eps;
Z3 = sin(r)./r;

figure(1);

subplot(2,2,1), contour(Z1)
title('sin(x) cos(y)')

subplot(2,2,2), contour(x,y,Z3)
title('sin(r)/r')

subplot(2,2,3), contour3(Z2,15)
title('x-0.5x^3+0.2y^2+1')

subplot(2,2,4), contour3(x,y,Z3)
title('sin(r)/r')

subplot(2,2,3); rotate3d;


figure(2);

subplot(2,2,1), mesh(Z1)
title('sin(x) sin(y)')

subplot(2,2,2), meshz(Z2)
title('x-0.5x^3+0.2y^2+1')

subplot(2,2,3), waterfall(Z2)
title('x-0.5x^3+0.2y^2+1')

subplot(2,2,4), meshc(Z3)
title('sin(r)/r')

%%%%%%%%%%%%%%%% end ex10.m %%%%%%%%%%%%%%%%%





