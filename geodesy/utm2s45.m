function [y,x] = utm2s45(N,E)
%UTM2S45  Conversion of UTM (N,E), zone 33 to (y,x) in system 45.
%         See Borre (1993): Landmaaling, Table 5.7

%Kai Borre 10-18-1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date:1998/10/23  $

y0 = 6100000;
x0 = 500000;
a = [ 4.2579934018E4;
      1.0004053319;
     -7.7039294168E-11;
      4.7460391692E-15;
      1.8259975589E-18;
     -2.0727959960E-22;
      7.0257273976E-27;
     -8.0529957628E-32];
b = [-4.2993547686E4;
     -1.4643915438E-3;
     -4.7898338026E-10;
      3.7432077762E-14;
     -4.8899686702E-18;
      2.0540300545E-22;
     -2.9421011013E-27;
      1.5387577766E-33];
ystar = N-y0; % shift of false origin
xstar = E-x0; 
p = 0;
q = 0;
for t = 8:-1:1
   z = p*ystar-q*xstar+a(t); % the real part is summed
   q = p*xstar+q*ystar+b(t); % the imaginary part is summed
   p = z;
end;
fprintf('\n y = %8.3f m and x = %8.3f m\n',p,-q);
%%%%%%%%%%%%%%%%%%%%%%% end utm2s45.m  %%%%%%%%%%%%%%%%%%%%%
