%GPS Lecture 4
%Date 12-Nov-1999
%
%CHECK_T  repairs over- and underflow of GPS time
%
%DMS2RAD  Conversion of degrees, minutes, and seconds to radians
%
%GET_EPH  The ephemerides contained in ephemeridesfile
%	       are reshaped into a matrix with 21 rows and
%	       as many columns as there are ephemerides.
%	            Typical call eph = get_eph('rinex_n.dat')
%
%SATPOS   Calculation of X,Y,Z coordinates at time t for given 
%         ephemeris eph
%
%SAV	    Plots of a satellite's acceleration and velocity
%
%SIGMA_TR Application of Variance propagation. Data from Example 9.2 
%         in Kai Borre (1995): GPS i landmaalingen
%
%WGS2ED50 Conversion of (phi,lambda) in WGS 84 to (Phi,Lambda) in ED 50
%
%%%%%%%%%%%%%%%% end contents.m  %%%%%%%%%%%%%%%%%%%%%%%%%

