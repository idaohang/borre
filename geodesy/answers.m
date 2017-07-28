%Exercises for Elementary Geodesy
%Written by Kai Borre
%October 26, 1998; revised October 10,1999

%Ex 1
%Compute the length of Equator for a=6378137 m.
a = 6378137.0;
equator = 2*pi*a

%Ex 2
%Consider a levelling line observed from point A at latitude
%56 degrees to point B at latitude 57 degrees at elevation 0
%meters. The points A and B are in the geoid; hence the
%potential difference is 0.
%Now increase the potential at latitude 56 degrees by C=1000 m^2/s^2;
%and we arrive at point E. We do a similar thing at latitude 57
%degrees and arrive at point D. The potential difference along
%the closed curve A-B-D-E is zero. When moving around this curve no
%work is performed: the field is conservative.
%However, if we convert the potential difference of 1000 m^2/s^2
%to meters at the two different latitudes we get different
%orthometric heights H_N!

%The WGS document mentioned in the following is
%World Geodetic System 1984, third edition 4 July 1997
%National Imagery and Mapping Agency. Technical report DMATR83502WGS84

%We assume WGS84 values
phi = [56; 57]*pi/180;	  % unit rad
%WGS84 document eq. (4-1)
gamma_e = 9.7803253359;
k = 6356752.3142*9.8321849378/(6378137*gamma_e)-1;
gamma = gamma_e*(1+k*(sin(phi)).^2)./  ...
		   sqrt(1-0.00669437999014.*(sin(2*phi)).^2);
f = 1/298.257223563;
%WGS84 document eq. (3-6) and (3-3)
m = (7292115*1.e-11)^2*6378137^3/(3986004.418*1.e8);
C = 1000;		  %~100 m
H_N = C./gamma; 
%WGS84 document eq. (4-3)
gamma_h = gamma.*(1-2*(1+f+m-2*f.*(sin(phi)).^2).*H_N/a+3*H_N.^2/a^2);
H_N = C./gamma_h;
fprintf('\n %12.5f\n',H_N(1), H_N(2));

%Ex 3
%The datum ED 50 is based on the international ellipsoid of 1924
%with a = 6378388 m. The datum WGS84 has a = 6378137 m.
%Adding the 251 m to the semi-axis a, changes the position
%phi = 57 degrees North, and lambda = 10 degrees East. But how?

%Lambda remains unchanges because of rotational symmetry, while
%the latitude changes an amount comparable to the 251 m.
%Later in the couse we learn how to compute the exact change.

%For the impatient student we quote the exact answer here
[X,Y,Z] = frgeod(6378388,297,57,10,0);
[phi,lambda,h] = togeod(6378137,298.257223563,X,Y,Z);
rad2dms((57-phi)*pi/180);
rad2dms((10-lambda)*pi/180);

%Ex 4
%The different position of ED50 with respect to WGS84 is decribed
%by the translational vector t=(-87,-98,-121). Can you describe the
%influence of t on the position phi = 57 degrees North, and
%lambda = 10 degrees East.

%Again for the impatient student we compute the exact answer
[pos(1,1),pos(2,1),pos(3,1)] = frgeod(6378137,298.257223563,57,10,0);
pos = pos+[-87;-98;-121];
[phi,lambda,h] = togeod(6378137,298.257223563,pos(1),pos(2),pos(3));
rad2dms((57-phi)*pi/180);
rad2dms((10-lambda)*pi/180);

%Ex 5
%Computing length of the meridian arc from Equator to latitude phi
dtr = pi/180;
phi = 56;
phi = phi*dtr;
%International Ellipsoid 1924
a = 6378388;
f = 1/297;
ex2 = (2-f)*f/(1-f)^2;
B = a*sqrt(1+ex2)*( ...
     (1-3/4*ex2+45/64*ex2^2-175/256*ex2^3+395/587*ex2^4)*phi ...
    +(-3/8*ex2+15/32*ex2^2-525/1024*ex2^3+316/587*ex2^4)*sin(2*phi)...
	     +(15/256*ex2^2-105/1024*ex2^3+79/587*ex2^4)*sin(4*phi)...
			+(-35/3072*ex2^3+105/4096*ex2^4)*sin(6*phi)...
				       +(59/24550*ex2^4)*sin(8*phi))

%Alternative procedure using the n-parameter.
%This yields a faster convergence. Implementation according to
%Koenig & Weise
n = f/(2-f);
Qm = a*(1+n^2/4+n^4/64)/(1+n);
B = Qm*(phi+(-3/2*n+9/16*n^3)*sin(2*phi)...
	  +(15/16*n^2-15/32*n^4)*sin(4*phi)...
	  -35/48*n^3*sin(6*phi)+(315/512*n^4)*sin(8*phi))

%Ex 6
phi = 1/3600;
phi = phi*dtr;
B = mer_quad(Qm,n,phi)	% 1'' at Equator
%curvature in prime vertical at Equator
N = a*sqrt(1+ex2)/sqrt(1+ex2*(cos(0))^2);
L = N*cos(0)*dtr/3600

%1'' at 57 degrees
B = mer_quad(Qm,n,(57+1/3600)*dtr)-mer_quad(Qm,n,57*dtr)
%curvature in prime vertical at 57 degrees
N = a*sqrt(1+ex2)/sqrt(1+ex2*(cos(57*dtr))^2);
L = N*cos(57*dtr)*dtr/3600

%Ex 7
%Direct problem. Given a point with phi = 57, lambda = 10
s = [1000; 10000; 30000];  % 1 km, 10 km, and 30 km.
b_acc = [];
l_acc = [];
az_acc = [];
for i = 1:3
   [b,l,az] = gauss_di([57 0 0],[10 0 0],[45 0 0],s(i),5);
   b_acc = [b_acc; b];
   l_acc = [l_acc; l];
   az_acc = [az_acc; az];
end

%Ex 8
%Inverse problem. Given the results form Ex 1: b,l,az
for i = 1:3
   [s,azi1,azi2] = gauss_in([57 0 0],[10 0 0],b_acc(i,:),l_acc(i,:),5);
end

%Ex 9
%Given a point with (phi,lambda,h) in WGS 84. Find the Cartesian
%coordinates (X,Y,Z)
[X,Y,Z]= geo2cart([57 0 0],[10 0 0],60,5);

%Ex 10
%Given a point with the (X,Y,Z) from Ex 3. Find the geographical
%coordinates of the point in WGS 84.

[phi,lambda,h] = cart2geo(X,Y,Z,5);

%Ex 11
%Given the same (X,Y,Z) as in Ex 4, but compute the geographical
%coordinates now relative to the international ellipsoid 1924.
[phi,lambda,h] = cart2geo(X,Y,Z,1);

%Comment the result!

%Ex 12
datumch

%Ex 13
datumch(6378388,56,10,0,0,0,-87,-98,-121,-251,1/297-1/298.257223563,0)

%Ex 14
[N,E] = cart2utm(3429122.662,604646.845,5325950.420,32);

%Ex 15
[phi,lambda] = utm2geo(6318036.28,560828.13,32);

%Ex 16
[N,E] = geo2utm(57,10,32);

%Ex 17
[phi,lambda] = utm2geo(6317972.081,560749.622,32);
%%%%%%%%%%%%%%%%%%%%%% end answers.m  %%%%%%%%%%%%%%%%%%%%%%%%%%
