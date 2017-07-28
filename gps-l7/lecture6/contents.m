%GPS Lecture 6
%Date 12-Nov-1999
%
%B_POINT  Prepares input to the Bancroft algorithm for finding
%	       a preliminary position of a receiver. The input is
%	       four or more pseudoranges and the coordinates of the
%	       satellites.
%
%BANCROFT Calculation of preliminary coordinates
%	       for a GPS receiver based on pseudoranges
%	       to 4 or more satellites. The ECEF
%	       coordinates (see function e_r_corr)
%	       are the first three elements of
%	       each row of B. The fourth element of each
%	       row of B contains the observed pseudorange.
%	       Each row pertains to one satellite.
%	       The pseudorange in the first row of B is
%	       used to descriminate between the two
%	       possible solutions.
%
%CHECK_T  repairs over- and underflow of GPS time
%
%E_R_CORR Returns rotated satellite ECEF coordinates
%	       due to Earth rotation during signal travel time
%
%FIND_EPH Finds the proper column in ephemeris array
%
%GET_EPH  The ephemerides contained in ephemeridesfile
%     	 are reshaped into a matrix with 21 rows and
%	       as many columns as there are ephemerides.
%	           Typical call eph = get_eph('rinex_n.dat')
%
%GET_RHO  Calculation of distance in ECEF system between
%	       satellite and receiver at time tR_RAW given the
%	       the pertinent ephemeris Eph.
%
%LORENTZ  Calculates the Lorentz inner product of the two
%         4 by 1 vectors x and y
%
%MAKEBASE We assume that you already have run the M-file bdata
%     	 with b-files from master and rover. The call would be
%	             bdata('masterfile','roverfile') or specifically
%		          bdata('b0810a94.076','b0005a94.076')
%	       That creates the file bdata.dat.
%	       Likewise we assume that the calls
%	             edata('masterfile','roverfile') and
%	             sdata('masterfile','roverfile')
%	       have created edata.dat and sdata.dat.
%	       Making Double Differences of Code and Phase Observations.
%	       Estimates ambiguities on L1 and L2, and baseline components.
%
%SATPOS   Calculation of X,Y,Z coordinates at time t for given 
%         ephemeris eph

%TOGEOD   Subroutine to calculate geodetic coordinates
%	       latitude, longitude, height given Cartesian
%	       coordinates X,Y,Z, and reference ellipsoid
%	       values semi-major axis (a) and the inverse
%	       of flattening (finv).
%         The units of linear parameters X,Y,Z,a must all agree 
%         (m,km,mi,ft,..etc). The output units of angular quantities
%         will be in decimal degrees (15.5 degrees not 15 deg 30 min).
%         The output units of h will be the same as the units of
%         X,Y,Z,a.
%
%TOPOCENT Transformation of vector dx into topocentric coordinate
%	       system with origin at X. Both parameters are 3 by 1 vectors.
%	       Output: D	vector length in units like the input
%		           Az	azimuth from north positive clockwise, degrees
%		           El	elevation angle, degrees
%
%TROPO	Calculation of tropospheric correction. The range correction
%        ddr in m is to be subtracted from pseudo-ranges and carrier
%        phases
%              sinel    sin of elevation angle of satellite
%              hsta	   height of station in km
%              p	      atmospheric pressure in mb at height hp
%              tkel	   surface temperature in degrees Kelvin at 
%                       height htkel
%              hum	   humidity in % at height hhum
%              hp	      height of pressure measurement in km
%              htkel    height of temperature measurement in km
%              hhum	   height of humidity measurement in km
%
%%%%%%%%% end contents.m  %%%%%%%%%%%%%%%%%%%
