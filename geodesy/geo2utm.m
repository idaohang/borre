function [NN,E] = geo2utm(B,L,zonenr)
%GEO2UTM  Conversion of geographical coordinates (B,L)
%         in degrees to UTM coordinates (N,E) in meters
%         in zone 'zonenr' (often = 32).

%Kai Borre 10-18-1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date:1998/10/23  $

% This impementation is based on O. Andersson & K. Poder (1981):
% Koordinattransformationer ved Geod\ae{}tisk Institut.
% Landinspekt\o{}eren, Vol. 30: 552--571 and Vol. 31: 76
% A further excellent reference (KW) is
% R. K\"onig & K.H. Weise (1951): Mathematische Grundlagen der
% h\"oheren Geod\"asie und Kartographie. Erster Band, Springer Verlag. *){$N+}

%Explanation of variables used:
%fl	            flattening of ellipsoid
%a	               semi major axis in m
%m0	            1-scale at central meridian; for UTM 0.0004
%Q_n              normalized meridian quadrant
%E0	            Easting of central meridian
%L0	            Longitude of central meridian
%bg(1) to bg(4)   constants for ellipsoidal geogr. to spherical geogr.
%gb(1) to gb(4)   constants for spherical geogr. to ellipsoidal geogr.
%gtu(1) to gtu(4) constants for ellipsoidal N, E to spherical N, E
%utg(1) to utg(4) constants for spherical N, E to ellipoidal N, E
%tolutm           tolerance for utm, 1.2E-10*meridian quadrant
%tolgeo           tolerance for geographical, 0.00040 second of arc 

fl = 1/297;
a = 6378388;
m0 = 0.0004;

% Normalized meridian quadrant, KW p. 50 (96), p. 19 (38b), p. 5 (21) 
n = fl/(2-fl);
m = n*n*(1/4+n*n/64);
w = (a*(-n-m0+m*(1-m0)))/(1+n);
Q_n = a+w;

% Easting and longitude of central meridian 
E0 = 500000;
L0 = (zonenr-30)*6-3;

% Check tolerance for reverse transformation 
tolutm = pi/2*1.2E-10*Q_n;
tolgeo = 0.000040;

% Coefficients of trigonometric series

%   ellipsoidal to spherical geographical , KW p. 186--187, (51)-(52)
%   bg(1) = n*(-2	+ n*(2/3    + n*(4/3	  + n*(-82/45))));
%   bg(2) = n^2*(5/3	+ n*(-16/15 + n*(-13/9)));
%   bg(3) = n^3*(-26/15 + n*34/21);
%   bg(4) = n^4*1237/630;
%   
%   spherical to ellipsoidal geographical , KW p. 190--191, (61)-(62)
%   gb(1) = n*(2        + n*(-2/3    + n*(-2	 + n*116/45)));
%   gb(2) = n^2*(7/3	+ n*(-8/5 + n*(-227/45)));
%   gb(3) = n^3*(56/15 + n*(-136/35));
%   gb(4) = n^4*4279/630;
%   
%   spherical to ellipsoidal N, E , KW p. 196, (69)
%   gtu(1) = n*(1/2	  + n*(-2/3    + n*(5/16     + n*41/180)));
%   gtu(2) = n^2*(13/48    + n*(-3/5 + n*557/1440));
%   gtu(3) = n^3*(61/240   + n*(-103/140));
%   gtu(4) = n^4*49561/161280;
%   
%   ellipsoidal to spherical N, E , KW p. 194, (65)
%   utg(1) = n*(-1/2	   + n*(2/3    + n*(-37/96	+ n*1/360)));
%   utg(2) = n^2*(-1/48    + n*(-1/15 + n*437/1440));
%   utg(3) = n^3*(-17/480 + n*37/840);
%   utg(4) = n^4*(-4397/161280);

% With fl:=1/297 we get 

bg(1) = -3.37077907E-3;
bg(2) = 4.73444769E-6;
bg(3) = -8.29914570E-9;
bg(4) = 1.58785330E-11;

gb(1) = 3.37077588E-3;
gb(2) = 6.62769080E-6;
gb(3) = 1.78718601E-8;
gb(4) = 5.49266312E-11;

gtu(1) = 8.41275991E-4;
gtu(2) = 7.67306686E-7;
gtu(3) = 1.21291230E-9;
gtu(4) = 2.48508228E-12;

utg(1) = -8.41276339E-4;
utg(2) = -5.95619298E-8;
utg(3) = -1.69485209E-10;
utg(4) = -2.20473896E-13;

% B,L refer to latitude and longitude. Southern latitude is negative 
% Ellipsoidal latitude, longitude to spherical latitude, longitude 

neg_geo = 'false';
if B < 0,
   neg_geo = 'true '
end;
Bg_r = abs(B*pi)/180;   %degrees to radians
res_clensin = clen_sin(bg,4,2*Bg_r);
Bg_r = Bg_r+res_clensin;
Lg_r = (L-L0)*pi/180;

% Spherical latitude, longitude to complementary spherical latitude
%      i.e. spherical N, E	    
cos_BN = cos(Bg_r);
Np = arg(cos(Lg_r)*cos_BN,sin(Bg_r));
Ep = arc_tanh(sin(Lg_r)*cos_BN);

% Spherical normalized N,E to ellipsoidal N, E 
Np = 2*Np;
Ep = 2*Ep;
[dN,dE] = clen_k_sin(gtu,4,Np,Ep);
Np = Np/2;
Ep = Ep/2;
Np = Np+dN;
Ep = Ep+dE;
NN = Q_n*Np;
E = Q_n*Ep+E0;
if neg_geo == 'true '
   NN = -N+20000000; 
end;
fprintf('\n Geographical coordinates transformed to UTM\n\n');
fprintf('phi_geo = %10.8f and lambda_geo = %11.8f\n',B,L);
fprintf('\n N  = %12.3f m and E  = %12.3f m\n',NN,E);

%-------------------------------------------------  

function result = clen_sin(ar,degree,argument);
cos_arg = 2*cos(argument);
hr1 = 0;
hr = 0;
for t = degree:-1:1
   hr2 = hr1;
   hr1 = hr;
   hr = ar(t)+cos_arg*hr1-hr2;
end;
result = hr*sin(argument);

%----------------------------------------------------

function [re,im] = clen_k_sin(ar,degree,arg_real,arg_imag)
sin_arg_r = sin(arg_real);
cos_arg_r = cos(arg_real);
sinh_arg_i = sinh(arg_imag);
cosh_arg_i = cosh(arg_imag);
r = 2*cos_arg_r*cosh_arg_i;
i = -2*sin_arg_r*sinh_arg_i;
hr1 = 0; hr = 0; hi1 = 0; hi = 0;
for t = degree:-1:1
   hr2 = hr1;
   hr1 = hr;
   hi2 = hi1;
   hi1 = hi;
   z = ar(t)+r*hr1-i*hi-hr2;
   hi =     i*hr1+r*hi1-hi2;
   hr = z;
end;
r = sin_arg_r*cosh_arg_i;
i = cos_arg_r*sinh_arg_i;
re = r*hr-i*hi;
im = r*hi+i*hr;

%----------------------------------------------------

function arc_tanhr = arc_tanh(x)
if abs(x) < 0.95
   arc_tanhr = log((1+x)/(1-x))/2.0;
else
   arc_tanhr = 1.0E30;
end;

%-----------------------------------------------------

function argr = arg(y,x)
argr = atan(x/y);
%%%%%%%%%%%%%%%%%%%%%%%%%% end geo2utm.m  %%%%%%%%%%%%%%%%%%%%%%
