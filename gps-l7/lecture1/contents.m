%GPS Lecture 1
%Version 1.0  11-Nov-1999
%
%B_CLOCK  Reading of binary P-code data as collected from Z-12 
%         receiver. Input of b-file from master. Typical call: 
%             bdata('b0810a94.076')
%
%CHECK_T  accounting for beginning or end of week crossover
%
%COMPTIME Reads receiver clock offset from a binary Ashtech 
%         observation file and plots it. See also b_clock
%
%DOY	    Calculation of day number of year. hour is split into hr,
%	       min, and sec
%
%GET_EPH  The ephemerides contained in ephemeridesfile are reshaped 
%         into a matrix with 21 rows and as many columns as there
%         are ephemerides. Typical call eph = get_eph('rinex_n.dat')
%
%GPS_TIME Conversion of Julian Day number to GPS week and Seconds 
%         of Week reckoned from Saturday midnight
%
%JULDAY   Conversion of date as given by
%	           y ... year (four digits)
%	           m ... month
%	           d ... day
%	           h ... hour and fraction hereof
%	       The conversion is only valid in the time span from March 
%         1900 to February 2100
%
%RECPOS   Least-squares searching for receiver position. Given 4 or 
%         more pseudoranges and ephemerides. Zoom on the plot to 
%         detect the search pattern! Idea to this script originates 
%         from Clyde C. Goad
%
%SATCONST Script for drawing GPS constellation in INERTIAL and 
%         ECEF frames
%
%SATPOS   Calculation of X,Y,Z coordinates at time t for given 
%         ephemeris eph
%
%SATPOSIN Calculation of X,Y,Z coordinates in an INERTIAL reference
%         frame at time t with given ephemeris eph
%
%%%%%%%%% end contents.m %%%%%%%%%
