%GPS Lecture 5
%Date 12-Nov-1999
%
%DD_COV   Computes the covariance matrix for double differenced
%         observations between r receivers and s satellites
%
%DE	    Analysis of the filter matrix of the four observation
%         filter N1 and N2. Transformation to the wide lane
%         ambiguity Nw and N1
%
%DMS2RAD  Conversion of degrees, minutes, and seconds to radians
%
%GPSVAR   Estimation of standard deviations for N1, N2, and N1-N2
%         when using a least-squares procedure in batch mode. There
%         are 2 code and 2 phase observations, and I is put equal
%         to zero.
%
%K_DD3	 Kalman Filter for Estimation of Ambiguities (with I = 0)
%	       Double differenced code and phase observations.
%	       SV is the satellite to be differenced with ref. sat. 26.
%	       The choices are: 2, 9, 16, 23, 27
%
%K_DD4    Kalman Filter for Estimation of Ambiguities. Double 
%         differenced code and phase observations
%
%ONE_WAY  Evaluation of one-way data. Observations from Z12 receiver
%         taken at master point -810 androver point -005 on March 17,
%         1994
%
%%%%%%%%%%%%%% contents.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
