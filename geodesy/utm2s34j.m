%function [y,x] = utm2s34j(N,E)
%UTM2S34J  Conversion of UTM (N,E), zone 32 to (y,x) in system 1934,
%          jysk-fynsk section. See Borre (1993): Landmaaling,
%		     Table 5.5

%Kai Borre 10-17-1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date:1998/10/23  $

y0 = 6231000;
x0 = 530000;
a = [ 1.9937114647E5;
      1.0001365021;
     -1.6486053260E-11;
     -1.8877421278E-17];
b = [-2.6538172775E5;
      1.9995520454E-2;
      3.4452862838E-10;
     -9.6740792784E-17];
y = input(' N = ');
x = input(' E = ');
ystar = y-y0; % shift to false origin
xstar = x-x0;
p = 0;
q = 0;
for t = 4:-1:1
   z = p*ystar-q*xstar+a(t);  % the real part is summed
   q = p*xstar+q*ystar+b(t);  % the imaginary part is summed
   p = z;
end;
fprintf('\n y = %10.3f m and x = %10.3f m\n',p,-q);
%%%%%%%%%%%%%% end utm2s34j.m  %%%%%%%%%%%%%%%%%%%
