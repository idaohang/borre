function [y,x] = utm2s34s(N,E)
%UTM2S45   Conversion of UTM (N,E), zone 32 to (y,x) in system 34,
%          section sjaelland. See Borre (1993): Landmaaling, Table 5.6

%Kai Borre 10-18-1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $ $$Date:1998/10/23  $

y0 = 6131000;
x0 = 678000;
a = [ 9.6432190004E4;
      9.9980074884E-1;
     -7.2077820180E-12;
      6.2510854612E-17];
b = [ -1.1935322229E5;
       1.9555849611E-2;
       2.1675209626E-9;
      -2.5215002720E-16];
ystar = N-y0; % shift to false origin
xstar = E-x0;
p = 0;
q = 0;
for t = 4:-1:1
   z = p*ystar-q*xstar+a(t);  % the real part is summed
   q = p*xstar+q*ystar+b(t);  % the imaginary part is summed
   p = z;
end;
fprintf('\n y =  %8.3f m and x = %8.3f m \n',p,-q);
%%%%%%%%%%%%%%%%%%%%%%% end utm2s34s.m	%%%%%%%%%%%%%%%

