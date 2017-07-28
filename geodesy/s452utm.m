function [p,q] = s452utm(y,x)
%S452UTM  Conversion of system 1945 to UTM, zone 33
%	      See Borre (1993): Landmaaling, Table 5.4

%Kai Borre 10-17-1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $ $Date:1998/10/23  $

y0 = 50000;
x0 = 50000;
a = [ 6.1074273382E6;
      9.9959968336E-1;
      1.2072957024E-10;
     -9.4593449880E-15;
      6.5994864072E-19;
      2.4410466467E-23;
     -2.8946629725E-27;
      8.0282195772E-32];
b  = [4.9300725696E5;
      1.4686762273E-3;
      5.4905973626E-10;
      2.1510833240E-15;
     -1.0163477904E-18;
      4.4726583308E-23;
     -1.1066966998E-27;
     -5.8658108998E-34];
ystar = y-y0;	     % shift to false origin
xstar = -(x-x0);
p = 0;
q = 0;
for t = 8:-1:1
   z = p*ystar-q*xstar+a(t);  % the real part is summed
   q = p*xstar+q*ystar+b(t);  % the imaginary part is summed
   p = z;
end;
fprintf('\n N =  %10.3f m and  E =  %10.3f m \n',p,q);
%%%%%%%%%%%%%%%%%%%%%%%%% end s452utm.m  %%%%%%%%%%%%%%%%%%%%%
