function [p,q] = s34s2utm(y,x)
%S34S2UTM Conversion of (y,x) in system 1934, sjaelland to (N,E) in
%         UTM, zone 32. See Kai Borre (1993): Landmaaling, Table 5.3

%Kai Borre 10-17-1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date:1998/10/23  $

y0 = 200000;
x0 = 200000;
a = [ 6.2329356122E6;
      9.9945508220E-1;
     -7.3897597436E-11;
     -5.2080308192E-17];
b = [ 5.9533534820E5;
     -1.9979251874E-2;
     -2.0717946709E-9;
      2.5703781294E-16];
ystar = y-y0;	    % shift to false origin
xstar = -(x-x0);
p = 0;
q = 0;
for t = 4:-1:1
   z = p*ystar-q*xstar+a(t); % the real part is summed
   q = p*xstar+q*ystar+b(t); % the imaginary part is summed
   p = z;
end;
fprintf('\n N =  %10.3f m and E = %10.3f m\n',p,q);
%%%%%%%%%%%%%%%%%%%%%%%%%% end s34s2utm.m  %%%%%%%%%%%%%%%%%%%%%%%%
