%ADJUSTMENT TOOLBOX
%Version 1.0 10-Nov-1999
%
%CLOWN	 Image compression by means of SVD.
%	       We illustrate low-rank approximations of a clown
%
%DD_COV   Computes the covariance matrix for double differenced
%         observations between r receivers and s satellites
%
%ELLCONF  Computes the eigenvalues and -vectors of a 2 x 2 covariance
%   	    matrix, prints the values of the semi major, semi minor
%	       axes and  the rotation angle of the pertinent confidence
%  	    ellipse and plots it.
%
%EX31	    Solution to Eksempel 3.1 in "Mindste kvadraters princip",
%  	    pages 63--65
%
%EX31FREE Solution to Eksempel 3.1, but free-type network
%
%EX68	    Solution to Eksempel 6.8
%
%EX69	    Solution to Eksempel 6.9
%
%FIXING   Filter version of ex31free. Demonstrates the impact of
%	       introducing constraints as observations with zero variance
%
%K_UPDATE Kalman update, one measurement per call.  Observation
%	       covariance R
%
%L1	    Script for solving Eksempel 1.22, 1.23, 1.24, and 1.25
%
%L2	    Script for computing "En rektangulær grund"
%
%L3	    Script for computation of partial derivatives of a
%         distance observation
%
%LEV	    Least squares estimation of heights in a levelling network
%	       as described by the following 3 data files:
%		    levfix.dat, contains row-wise 2 columns
%				point#	 elevation
%			 levfree.dat, contains row-wise 1 column
%				 point#
%			 levobs.dat, contains row-wise 4 columns
%			    from-#    to-#	  HDIFF      std. dev. of HDIFF
%	       All units in meters
%
%OPG65	 Solution to Opgave 6.5
%
%OPG67	 Solution to Opgave 6.7
%
%RLS	    Recursive Least Squares. A is the coefficient matrix,
%	       b the observations, and Sigma a vector containing the
%	       diagonal entries of the covariance matrix for the problem.
%
%SET	    All direction observations are entered into a matrix.
%	       A column contains the direction as observed in each set.
%	       The zero direction is omitted.
%	       Output is the estimated direction, the orientation unknown
%	       and the standard deviation of one direction as observed by
%	       one set.
%
%SUR_DD_C Surface plot of the covariance matrix for double differenced
%	       observations
%
%%%%%%%%%%%%%%%%%%% end contents.m  %%%%%%%%%%%%%%%%%%%
