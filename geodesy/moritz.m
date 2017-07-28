function [A1,A2,s] = moritz(phi1d,phi1m,phi1s,lambda1d,lambda1m,...
                            lambda1s,phi2d,phi2m,phi2s,...
                            lambda2d,lambda2m,lambda2s)
%MOTITZ  Solution of the Reverse Geodetic Problem as given in
%	      Moritz, Helmut (1959): Eine direkte L\"osung der
%	      zweiten Hauptaufgabe auf dem Rotationsellipsoid
%	      f\"ur beliebige Entfernungen.
%	      Zeitschrift f\"ur Vemessungswesen, 84: 453--457

%Kai Borre, January 1, 1999
%Copyright (c) by Kai Borre
%$Revision 1.0 $  $Date:1999/01/01  $

% Ellipsoidal parameters for WGS84
a = 6378137;	         	% semi-major axis, meters
f = 1/298.257223563;	      % flattening
ex2 = (2-f)*f/(1-f)^2;	   % second eccentricity squared
c = sqrt(1+ex2)*a;	      % polar radius
dtr = pi/180;              % degrees to radians

if nargin == 0
   phi1  = 55.66167444*dtr;
   lambda1 = 0.300938611*dtr;
   phi2 = 55.81532861*dtr;
   lambda2 = 0.41109333*dtr;
end 
phi1 = dms2rad(phi1d,phi1m,phi1s);
lambda1 = dms2rad(lambda1d,lambda1m,lambda1s);
phi2 = dms2rad(phi2d,phi2m,phi2s);
lambda2 = dms2rad(lambda2d,lambda2m,lambda2s);
l = lambda2-lambda1;

% Computation of C0 according to (7)
ar =  tan(phi1)^2 + tan(phi2)^2 - 2*tan(phi1)*tan(phi2)*cos(l);
if l ~= 0
   C0 = atan(sqrt(ar)/sin(l));
else
   C0 = 0;
end

%first definitions of auxiliary quantities
u1 = asin(sin(phi1)/sin(C0));
u2 = asin(sin(phi2)/sin(C0));
du = u2-u1;
dtan1u = tan(u2)-tan(u1);
dtan3u = tan(u2)^3-tan(u1)^3;
dtan5u = tan(u2)^5-tan(u1)^5;
dtan7u = tan(u2)^7-tan(u1)^7;

dsin2u = sin(2*u2)-sin(2*u1);
dsin4u = sin(4*u2)-sin(4*u1);
dsin6u = sin(6*u2)-sin(6*u1);

% definition of coefficients in (9)
a(1,1) = -dtan1u/sin(C0);
a(1,2) = cot(C0)*(2*dtan1u+dtan3u)/(2*sin(C0));
a(1,3) = cot(C0)^2*((-6-2*tan(C0)^2)*dtan1u+ ...
          (-7-tan(C0)^2)*dtan3u-3*dtan5u)/(6*sin(C0));
a(1,4) = cot(C0)^3*((24+16*tan(C0)^2)*dtan1u+ ...
          (48+20*tan(C0)^2)*dtan3u+...
          (45+9*tan(C0)^2)*dtan5u+15*dtan7u)/(24*sin(C0));

a(2,1) = cos(C0)*cot(C0)*(tan(C0)^2*du+dtan1u)/2;
a(2,2) = cos(C0)*cot(C0)^2*(tan(C0)^2*du+ ...
          (-2-3*tan(C0)^2)*dtan1u-dtan3u)/4;
a(2,3) = cos(C0)*cot(C0)^3*(-tan(C0)^4*du+...
          (6+8*tan(C0)^2+3*tan(C0)^4)*dtan1u+ ...
          (7+6*tan(C0)^2)*dtan3u+3*dtan5u)/12;
a(3,1) = 3*cos(C0)*cot(C0)*((-tan(C0)^2-3*sin(C0)^2)*du/2-...
          cos(C0)^2*dtan1u-tan(C0)^2*sin(C0)^2*dsin2u/4)/8;
a(3,2) = 3*cos(C0)*cot(C0)^2*((5*tan(C0)^2-9*sin(C0)^2)*du/2+...
          (6-4*cos(C0)^2)*dtan1u+cos(C0)^2*dtan3u-...
          tan(C0)^2*sin(C0)^2*dsin2u/4)/16;
a(4,1) = 5*cos(C0)*cot(C0)*(tan(C0)^2*(24-36*sin(C0)^2+...
          15*sin(C0)^4)*du/8+cos(C0)^4*dtan1u+...
          tan(C0)^2*(3*sin(C0)^2-2*sin(C0)^4)*dsin2u/4+...
          tan(C0)^2*sin(C0)^4*dsin4u/32)/16;
dl1 = -cos(C0)*du/2;
dl2 = 3*cos(C0)*((2-sin(C0)^2)*du+sin(C0)^2*dsin2u/2)/16;
dl3 = 5*cos(C0)*((-8+8*sin(C0)^2-3*sin(C0)^4)*du/2+...
          (-2*sin(C0)^2+sin(C0)^4)*dsin2u-sin(C0)^4*dsin4u/8)/64;
dl4 = 35*cos(C0)*((16-24*sin(C0)^2+18*sin(C0)^4-...
          5*sin(C0)^6)*du+(48*sin(C0)^2-48*sin(C0)^4+...
          15*sin(C0)^6)*dsin2u/4+(6*sin(C0)^4-3*sin(C0)^6)*...
          dsin4u/4+sin(C0)^6*dsin6u/12)/2048; 

% Solution of system (8)         
C1 = -dl1/a(1,1);
C2 = -(a(1,2)*C1^2+a(2,1)*C1+dl2)/a(1,1);
C3 = -(a(1,2)*2*C1*C2+a(1,3)*C1^3+a(2,1)*C2+...
          a(2,2)*C1^2+a(3,1)*C1+dl3)/a(1,1);
C4 = -(a(1,2)*(2*C1*C3+C2^2)+a(1,3)*3*C1^2*C2+a(1,4)*C1^4+...
          a(2,1)*C3+a(2,2)*2*C1*C2+a(2,3)*C1^3+a(3,1)*C2+...
          a(3,2)*C1^2+a(4,1)*C1+dl4)/a(1,1);
C = C0+ex2*C1+ex2^2*C2+ex2^3*C3+ex2^4*C4;

% Computation of new variables
U1 = asin(sin(phi1)/sin(C));
U2 = asin(sin(phi2)/sin(C));

dU = U2-U1;
dsin2U = sin(2*U2)-sin(2*U1);
dsin4U = sin(4*U2)-sin(4*U1);
dsin6U = sin(6*U2)-sin(6*U1);
dsin8U = sin(8*U2)-sin(8*U1);

% Computation of azimuths
alpha1 = atan(cot(C)*sqrt(1+ex2*cos(phi1)^2)/cos(U1));
alpha2 = atan(cot(C)*sqrt(1+ex2*cos(phi2)^2)/cos(U2));
A1 = rad2dms(alpha1);
A2 = rad2dms(alpha2);

% Computation of coefficients according to (12)
s0 = c*dU; 
s1 = c*((-4+sin(C)^2)*dU/4-3*sin(C)^2*dsin2U/8);
s2 = c*((64-32*sin(C)^2+13*sin(C)^4)*dU/64+...
         (24*sin(C)^2-9*sin(C)^4)*dsin2U/32+15*sin(C)^4*dsin4U/256);
s3 = c*((-256+192*sin(C)^2-156*sin(C)^4+45*sin(C)^6)*dU/256+ ...
         (-1152*sin(C)^2+864*sin(C)^4-237*sin(C)^6)*dsin2U/1024+...         
         (-180*sin(C)^4+75*sin(C)^6)*dsin4U/1024-...
         35*sin(C)^6*dsin6U/3072);      
s4 = c*((16384-16384*sin(C)^2+19968*sin(C)^4-11520*sin(C)^6+...
         2577*sin(C)^8)*dU/16384+(6144*sin(C)^2-6912*sin(C)^4+...
         3792*sin(C)^6-819*sin(C)^8)*dsin2U/4096+...
         (5760*sin(C)^4-4800*sin(C)^6+1245*sin(C)^8)*dsin4U/16384+...
         (560*sin(C)^6-245*sin(C)^8)*dsin6U/12288+...
         315*sin(C)^8*dsin8U/131072);
s = s0+ex2*s1+ex2^2*s2+ex2^3*s3+ex2^4*s4;

%--------------------------------------------

function result = dms2rad(deg,min,sec)
%DMS2RAD Conversion of degrees, minutes, and seconds to radians

arg = deg+min/60+sec/3600;
result = arg*pi/180;

%-------------------------------

function result = rad2dms(arg)
%RAD2DMS Conversion of radians to degrees, minutes, and seconds

arg = arg*180/pi;
result = zeros(1,3);
result(1) = fix(arg);
if result(1) == 0
   result(2) = fix(arg*60);
else
   result(2) = fix(rem(arg,result(1))*60);
end
result(3) = (arg-result(1)-result(2)/60)*3600;
%fprintf('\n %3.0f %3.0f %10.6f\n',result(1),result(2),result(3))

%%%%%%%%%%%%%%%%% end moritz.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%
