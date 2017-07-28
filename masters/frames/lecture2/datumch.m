function change = datumch(a,phi,lambda,ex,ey,ez,...
                  					 tx,ty,tz,da,df,dk)
%DATUMCH Computes the change of latitude (degrees),
%	 longitude (degrees) and height (meters) by rotating
%	 the ellipsoid the angles ex,ey,ez (degrees), by
%	 translating the ellipsoid tx,ty,tz (meters), by
%	 changing the semi-major axis a (meters) by
%	 da (meters) and the flattening by df (dimensionless),
%	  and by changing the scale by dk (dimensionless).

%Kai Borre 10-16-1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date:1998/10/23  $

if nargin == 0
   a = 6378388; %Transformation from ED50 to WGS84
   phi = 57;	 %degrees
   lambda = 10; %degrees
   ex = 0.0;
   ey = 0.0;
   ez = 0.0;
   tx = -87.0;	 %meters
   ty = -98.0;	 %meters
   tz = -121.0;  %meters
   da = -251.0;  %meters
   df = 1/297-1/298.257223563;
   dk = 0.0;
end

dtr = pi/180;
phi = phi*dtr;
lambda = lambda*dtr;
rot = [ex; ey; ez]*dtr;
trans = [tx; ty; tz];
ellch = [da; df];

F = [-sin(lambda) -sin(phi)*cos(lambda) cos(phi)*cos(lambda);
      cos(lambda) -sin(phi)*sin(lambda) cos(phi)*sin(lambda);
	       0	             cos(phi)		            sin(phi)];

DR = [-a*sin(phi)*cos(lambda) -a*sin(phi)*sin(lambda) a*cos(phi);
	       a*sin(lambda)	      -a*cos(lambda)	                  0;
	          0		                   0		                    0];

FTG = [0		                   0;
       0 a*sin(phi)*(cos(phi))^3;
       1	       -a*(sin(phi))^2];

change = DR*rot-F'*trans-[0;0;a]*dk-FTG*ellch;
ch = change(1,1)/(a*cos(phi)); % we switch components 1 and 2
change(1,1) = change(2,1)/a;
change(2,1) = ch;
rad2dms(change(1,1));
rad2dms(change(2,1));
%fprintf('\n\n Change in latitude [dd mm ss.sssss]: %4.0f %4.0f %12.6f',...
%fprintf('\n Change in longitude [dd mm ss.sssss]: %4.0f %4.0f %12.6f',...
%fprintf('\n Change in height [m]: %10.4f\n', change(3,1));
%%%%%%%%%%%%%%%%%%%% end datumch.m  %%%%%%%%%%%%%%%%




