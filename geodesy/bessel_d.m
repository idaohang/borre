function bessel_d(phi1d,phi1m,phi1s,l1d,l1m,...
                                    l1s,A1d,A1m,A1s,s12,a,finv)
%BESSEL_D Solution of the direct geodetic problem according to
%         the Bessel-Helmert method as described in Bodem\"uller.
%         Given a point with coordinates (phi1, l1) and 
%         a gedesic with azimuth A1 and length s12 from here. The 
%         given reference ellipsoid has semi-major axis a and 
%         inverse flattening finv. The coordinates and the azimuth
%         have the format degree, minute, and second
%         with decimals. 
%
%  	    Bodem\"uller, H.(1954): Die geod\"atischen Linien des 
%	       Rotationsellipsoides und die L\"osung der geod\"atischen
%	       Hauptaufgaben f\"ur gro\ss{}e Strecken unter 
%         besonderer Ber\"ucksichtigung der Bessel-Helmertschen
%         L\"osungsmethode.
%	       Deutsche Geod\"atische Kommission, Reihe B, Nr. 13.

%Kai Borre, January 10, 1999
%Copyright (c) by Kai Borre
%$Revision 1.0 $  $Date:1999/01/10  $

dtr = pi/180;        % degrees to radians
if nargin == 0
   phi1  = 50*dtr;
   l1 =    10*dtr;
   A1 =   140*dtr;
   s12 = 15000000;   % m
   a = 6378388;
   finv = 297;
else
   phi1 = dms2rad(phi1d,phi1m,phi1s);
   l1 = dms2rad(l1d,l1m,l1s);
   A1 = dms2rad(A1d,A1m,A1s);
end

f = 1/finv; 
ex1 = (2-f)*f;             % first eccentricity squared
ex2 = (2-f)*f/(1-f)^2;	   % second eccentricity squared
%The implementation follows pages 38-39 in Bodem\"uller.
%reduced latitude beta1
tanbeta1 = (1-f)*tan(phi1); 
sigma11 = - atan(cos(A1)/tanbeta1);
beta1 = atan(tanbeta1); 
%largest reduced latitude for the geodesic
cosbetam = cos(beta1)*sin(A1);
if cosbetam == 0
   cosbetam = 1.e-12;
end
lambda1 = atan(tan(sigma11)/cosbetam);
sinbetam = sqrt(1-cosbetam^2);
K1 = (sqrt(1+ex2*sinbetam^2)-1)/(sqrt(1+ex2*sinbetam^2)+1); % (52)
s11 = sigma11+K1*sin(2*sigma11)/2-K1^2*sin(4*sigma11)/16;   % (61)
deltas1 = s12*(1-K1)/((1-f)*a*(1+K1^2/4)); % (58')
s1 = s11+deltas1/2;
deltasigma1 = deltas1-(K1-9*K1^3/16)*cos(2*s1)*sin(deltas1)+...
              5*K1^2*cos(4*s1)*sin(2*deltas1)/8-...
              29*K1^3*cos(6*s1)*sin(3*deltas1)/48;    % (62)
sigma12 = sigma11+deltasigma1;             % (63)
sinsigma12 = sin(sigma12);
tansigma12 = tan(sigma12);
lambda2 = atan(tansigma12/cosbetam);
beta2 = acos(sin(sigma12)/sin(lambda2));
A2 = -acos(tansigma12*tan(beta2));
if A2 < 0, A2 = A2+pi; end
phi2 = atan(tan(beta2)/(1-f));
l2 = l1 +(lambda2-lambda1) ...
       -(2-f)*f*cosbetam*((1+f/(2-f)-K1/2-K1^2/4)*deltasigma1 ...
       -K1*cos(2*s1)*sin(deltasigma1)/2 ...
       +K1^2*cos(4*s1)*sin(2*deltasigma1)/8)/2;    % (57)
if l2 < 0, l2 = l2+pi; end
disp('A2');
rad2dms(A2);        
disp('phi2')
rad2dms(phi2);        
disp('lambda2')
rad2dms(l2);

%----------------------------------------------

function result = dms2rad(deg,min,sec);
% Conversion of degrees, minutes, and seconds to radians

neg_arg = 'FALSE';
if deg < 0
   neg_arg = 'TRUE ';
   deg = -deg;
end
arg = deg+min/60+sec/3600;
result = arg*pi/180;
if neg_arg == 'TRUE ';
   result = -result;
end

%------------------------------------------

function result = rad2dms(arg)
%RAD2DMS Conversion of radians to degrees, minutes, and seconds

neg_arg = 'FALSE';
if arg < 0
   neg_arg = 'TRUE ';
   arg = -arg;
end

arg = arg*180/pi;
result = zeros(1,3);
result(1) = fix(arg);
if result(1) == 0
   result(2) = fix(arg*60);
else
   result(2) = fix(rem(arg,result(1))*60);
end
result(3) = (arg-result(1)-result(2)/60)*3600;
if neg_arg == 'TRUE '
   result(1) = -result(1);
end
fprintf('   %3.0f %2.0f %8.6f\n',result(1),result(2),result(3))

%%%%%%%%%%%%%%%%% end bessel_d.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%

