function bessel_i(phi1d,phi1m,phi1s,l1d,l1m,...
                      l1s,phi2d,phi2m,phi2s,l2d,l2m,l2s,a,finv)
%BESSEL_I Solution of the inverse geodetic problem according to
%         the Bessel-Helmert method as described in Bodem\"uller.
%         Given two points with coordinates (phi1, l1) and 
%         (phi2,l2). Per definition we always have l2 > l1.
%         The given reference ellipsoid has semi-major
%         axis a and inverse flattening finv. The coordinates
%         have the format degree, minute, and second with decimals. 
%
%  	    Bodem\"uller, H.(1954): Die geod\"atischen Linien des 
%	       Rotationsellipsoides und die L\"osung der geod\"atischen
%	       Hauptaufgaben f\"ur gro\ss{}e Strecken unter 
%         besonderer Ber\"ucksichtigung der Bessel-Helmertschen
%         L\"osungsmethode.
%	       Deutsche Geod\"atische Kommission, Reihe B, Nr. 13.

%Kai Borre, January 12, 1999
%Copyright (c) by Kai Borre
%$Revision 1.0 $  $Date:1999/01/12  $

if nargin == 0
   phi1  = dms2rad(50,0,0);
   l1 = dms2rad(10,0,0);
   phi2 = dms2rad(-62,57,3.203824);
   l2 = dms2rad(105,5,38.299430);
   a = 6378388;
   finv = 297;
else
   phi1 = dms2rad(phi1d,phi1m,phi1s);
   l1 = dms2rad(l1d,l1m,l1s);
   phi2 = dms2rad(phi2d,phi2m,phi2s);
   l2 = dms2rad(l2d,l2m,l2s);
end
dtr = pi/180;
f = 1/finv; 

%The implementation follows pages 41-42 in Bodem\"uller.
%Reduced latitudes
beta1 = atan((1-f)*tan(phi1)); 
beta2 = atan((1-f)*tan(phi2));
lambda_dif = (l2-l1)/2;

for i = 1:5 % iteration for betam
   lambda_sum = acot(-tan(lambda_dif)* ... 
                    sin(beta2+beta1)/sin(beta2-beta1)); % (66)
   if lambda_sum < 0
      lambda_sum = lambda_sum+pi;
   end
   lambda2 = lambda_sum+lambda_dif;
   lambda1 = lambda_sum-lambda_dif;
   betam = atan(tan(beta1)/cos(lambda1));            % (67)
   s11 = acos(sin(beta1)/sin(betam));
   s21 = acos(sin(beta2)/sin(betam));
   deltasigma1 = s21-s11;
   s1 = (s21+s11)/2;
   sqK = sqrt(1+(2-f)*f*(sin(betam))^2/((1-f)^2));
   K1 = (sqK-1)/(sqK+1);                             % (52)
   lambda_dif = (l2-l1)/2 ...
                   +(2-f)*f*cos(betam)*(...
                  (1+f/(2-f)-K1/2-K1^2/4)*deltasigma1 ...
                  -K1*cos(2*s1)*sin(deltasigma1)/2 ...
                  +K1^2*cos(4*s1)*sin(2*deltasigma1)/8)/4; % (57)                   
end

A1 = acot(-tan(betam)*sin(s11));
if A1 < 0
   A1 = A1+pi;
end
A2 = acot(-tan(betam)*sin(s21));
if A2 < 0
   A2 = A2+pi;
end
disp('A1');
rad2dms(A1);        
disp('A2')
rad2dms(A2);        
deltas1 = deltasigma1+(K1-3*K1^3/8)*cos(2*s1)*sin(deltasigma1) ...
                       -K1^2*cos(4*s1)*sin(2*deltasigma1)/8; % (69)
deltas = deltas1*(1-f)*a*(1+K1^2/4)/(1-K1);
fprintf('\ns   %8.4f m\n', deltas)

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

%%%%%%%%%%%%%%%%% end bessel_i.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%
