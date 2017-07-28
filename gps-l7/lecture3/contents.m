
%GPS Lecture 3
%Date 12-Nov-1999
%
%CHECK_T  repairs over- and underflow of GPS time
%
%FIND_EPH Finds the proper column in ephemeris array
%
%FRGEOD   Subroutine to calculate Cartesian coordinates X,Y,Z
%	       given geodetic coordinates latitude, longitude (east),
%	       and height above reference ellipsoid along with
%	       reference ellipsoid values semi-major axis (a) and
%	       the inverse of flattening (finv).
%	       The units of linear parameters h,a must agree 
%         (m,km,mi,..etc). The input units of angular quantities
%         must be in decimal degrees. The output units of X,Y,Z
%         will be the same as the units of h and a.
%
%GET_EPH  The ephemerides contained in ephemeridesfile are reshaped
%         into a matrix with 21 rows and as many columns as there 
%         are ephemerides.
%	       Typical call eph = get_eph('rinex_n.dat')
%
%RAD2DMS  Conversion of radians to degrees, minutes, and seconds
%
%RECPOS   Given 4 or more pseudoranges and ephemerides.
%	       Zoom on the plot to detect the search pattern!
%         Idea to this script originates from Clyde C. Goad
%
%SATPOS   Calculation of X,Y,Z coordinates at time t for given
%         ephemeris eph
%
%TOGEOD   Subroutine to calculate geodetic coordinates latitude,
%         longitude, height given Cartesian coordinates X,Y,Z,
%         and reference ellipsoid values semi-major axis (a)
%         and the inverse of flattening (finv).
%         The units of linear parameters X,Y,Z,a must all agree
%         (m,km,mi,ft,..etc). The output units of angular quantities
%         will be in decimal degrees (15.5 degrees not 15 deg 30 min).  
%         The output unit of h will be the same as the units of
%         X,Y,Z,a.
%
%%%%%%%% end contents.m  %%%%%%%%%%%%%%%%%%%%%%
