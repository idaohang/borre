function [N,E] = s34j2utm(y,x)
%S34J2UTM  Conversion of (y,x) in system 1934 to (N,E) in UTM,
%          zone 32. See Borre (1993): Landmaaling, Table 5.2

%Kai Borre 10-17-1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date:1998/10/23  $

y0 = 200000;
x0 = 200000;
a = [6.2329350512E6;
     9.9950870832E-1;
    -2.2732469700E-11;
     2.6264691305E-17];
b = [5.9533558288E5;
    -1.9984265541E-2;
    -3.3921905837E-10;
     9.4810550416E-17];
ystar = y-y0;	 % shift to false origin
xstar = -(x-x0);
p = 0;
q = 0;
for t = 4:-1:1,
   z = p*ystar-q*xstar+a(t);  % the real part is summed
   q = p*xstar+q*ystar+b(t);  % the imaginary part is summed
   p = z;
end;
fprintf('\n N = %10.3f m and E = %10.3f m\n',p,q);
%%%%%%%%%%%%%%%% end s34j2utm.m %%%%%%%%%%%%%%%%
