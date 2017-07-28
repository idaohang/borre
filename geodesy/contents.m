%Geodesy Toolbox
%
%Version 1.1 12-Jan-1999
%
%Directory: geodesy
%
%BESSEL_D Solution of the direct geodetic problem according to
%	       the Bessel-Helmert method as described in Bodem\"uller.
%	       Given a point with coordinates (phi1, l1) and
%	       a gedesic with azimuth A1 and length s12 from here. The
%	       given reference ellipsoid has semi-major axis a and
%	       inverse flattening finv. The coordinates and the azimuth
%	       have the format degree, minute, and second
%	       with decimals.
%
%	       Bodem\"uller, H.(1954): Die geod\"atischen Linien des
%	       Rotationsellipsoides und die L\"osung der geod\"atischen
%	       Hauptaufgaben f\"ur gro\ss{}e Strecken unter
%	       besonderer Ber\"ucksichtigung der Bessel-Helmertschen
%	       L\"osungsmethode.
%	       Deutsche Geod\"atische Kommission, Reihe B, Nr. 13.
%
%BESSEL_I Solution of the inverse geodetic problem according to
%     	 the Bessel-Helmert method as described in Bodem\"uller.
%         Given two points with coordinates (phi1, l1) and
%	       (phi2,l2). Per definition we always have l2 > l1.
%	       The given reference ellipsoid has semi-major
%	       axis a and inverse flattening finv. The coordinates
%	       have the format degree, minute, and second with decimals.
%
%CART2GEO Conversion of Cartesian coordinates (X,Y,Z) to geographical
%	       coordinates (phi,lambda,h) on a selected reference ellipsoid
%	       Choices i of Reference Ellipsoid
%				1. International Ellipsoid 1924
%			   2. International Ellipsoid 1967
%				3. World Geodetic System 1972
%				4. Geodetic Reference System 1980
%				5. World Geodetic System 1984
%
%CART2UTM Transformation of (X,Y,Z) to (N,E,U)in UTM, zone 'zone'
%
%DATUMCH  Computes the change of latitude (degrees),
%	       longitude (degrees) and height (meters) by rotating
%	       the ellipsoid the angles ex,ey,ez (degrees), by
%	       translating the ellipsoid tx,ty,tz (meters), by
%	       changing the semi-major axis a (meters) by
%	       da (meters) and the flattening by df (dimensionless),
%     	 and by changing the scale by dk (dimensionless).
%
%DMS2RAD  Conversion of degrees, minutes, and seconds to radians
%
%FRGEOD   Computing Cartesian coordinates X,Y,Z
%	       given geodetic coordinates latitude, longitude (east),
%	       and height above reference ellipsoid along with
%	       reference ellipsoid values semi-major axis (a)
%	       and the inverse of flattening (finv)
%	       The units of linear parameters h,a must agree
%	       (m,km,mi,..etc). The input units of angular quantities
%	       must be in decimal degree. The output units of X,Y,Z
%	       will be the same as the units of h and a.
%
%GAUSS_DI The geodetic direct problem is solved iteratively
%   		 by means of Gauss' mid-latitude formulas.
%	       Given the coordinates (phi1, lambda1) of a point,
%	       an azimuth az1 and distance s in meters.
%	       The format of phi1,lambda1 and az1 is [degrees minutes
%	       seconds]
%	       The choices i of Reference Ellipsoid are the following
%			      1. International Ellipsoid 1924
%			      2. International Ellipsoid 1967
%	            3. World Geodetic System 1972
%	 		      4. Geodetic Reference System 1980
%			      5. World Geodetic System 1984
%	       We look for the coordinates (phi2,lambda2) and
%	       the reverse azimuth az2.
%
%GAUSS_IN The geodetic inverse problem solved by means of
%	       Gauss' mid-latitude formulas.
%	       Given the coordinates (phi1, lambda1), and (phi2, lambda2)
%	       of two points in the format [degrees minutes seconds].
%	       The choices i of Reference Ellipsoid are the following:
%		    	  1. International Ellipsoid 1924
%			     2. International Ellipsoid 1967
%			     3. World Geodetic System 1972
%			     4. Geodetic Reference System 1980
%				  5. World Geodetic System 1984
%	       Unknowns are the distance s (in meters) and the mutual
%	       azimuths (in degrees).
%
%GEO2CART Conversion of geographical coordinates (phi,lambda,h)
%	       to Cartesian coordinates (X,Y,Z).
%	       Format for phi and lambda [degrees minutes seconds]
%	       h, X, Y, and Z are in meters
%
%GEO2UTM  Conversion of geographical coordinates (B,L)
%	       in degrees to UTM coordinates (N,E) in meters
%	       in zone 'zonenr' (often = 32).
%
%MER_QUAD Computes the normalized meridian from Equator to
%	       latitude phi
%
%RAD2DMS  Conversion of radians to degrees, minutes, and seconds
%
%S34J2UTM Conversion of (y,x) in system 1934 to (N,E) in UTM, zone 32. See
%	       Borre (1993): Landmaaling, Table 5.2
%
%S34S2UTM Conversion of (y,x) in system 1934, sjaelland to (N,E) in
%	       UTM, zone 32. See Kai Borre (1993): Landmaaling, Table 5.3
%
%S452UTM  Conversion of system 1945 to UTM, zone 33
%	       See Borre (1993): Landmaaling, Table 5.4
%
%TOGEOD   Computing geodetic coordinates
%	       latitude, longitude, height given Cartesian
%	       coordinates X,Y,Z, and reference ellipsoid
%	       values semi-major axis (a) and the inverse
%	       of flattening (finv).
%     	 The units of linear parameters X,Y,Z,a must all agree
%	       (m,km,mi,ft,..etc). The output units of angular quantities
%	       will be in decimal degrees (15.5 degrees not 15 deg 30 min).
%	       The output units of h will be the same as the units of
%	       X,Y,Z,a.
%
%UTM2GEO  Conversion of (N,E) (in meters) in UTM, zone 'zonenr' to
%	       geographical coordinates (phi,lambda) (in degrees)
%
%UTM2S34J Conversion of UTM (N,E), zone 32 to (y,x) in system 1934,
%         section jysk-fynsk. See Borre (1993): Landmaaling,
%	       Table 5.5
%
%UTM2S34S Conversion of UTM(N,E), zone 32 to (y,x) in system 1934,
%	       section sjaelland. See Borre (1993): Landmaaling,
%   	    Table 5.6
%
%UTM2S45  Conversion of UTM (N,E), zone 33 to (y,x) in system 1945.
%	       See Borre (1993): Landmaaling, Table 5.7
%
%WGS2ED50 Conversion of (phi,lambda) in WGS 84 to (Phi,Lambda)
%	       in ED 50
%%%%%%%%%%%%%%%%%%%%% end contents.m  %%%%%%%%%%%%%%%%%%%%%%%
