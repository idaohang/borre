%GPS Lecture 2
%Date 12-Nov-1999
%
%CHECK_T  accounting for beginning or end of week crossover
%
%EDATA	 Reads a binary ephemeris file and stores it in
%	       a matrix with 21 rows; column number is
%	       the number of ephemerides. eph is stored in efile.
%	       Typical call: edata('e0810a94.076','edata.dat')
%
%FIND_EPH Finds the proper column in ephemeris array
%
%GET_EPH  The ephemerides contained in ephemeridesfile
%	       are reshaped into a matrix with 21 rows and
%  	    as many columns as there are ephemerides.
%	       Typical call eph = get_eph('rinex_n.dat')
%
%GPS_TIME Conversion of Julian Day number to GPS week and
%	       Seconds of Week reckoned from Saturday midnight
%
%JULDAY   Conversion of date as given by
%	            y ... year (four digits)
%	            m ... month
%              d ... day
%	            h ... hour and fraction hereof
%	       The conversion is only valid in the time span
%	       from March 1900 to February 2100
%
%RUN	    Script for Lecture 2
%
%SATPOS   Calculation of X,Y,Z coordinates at time t for given
%         ephemeris eph
%
%%%%%%%%% end contents.m %%%%%%%%%
