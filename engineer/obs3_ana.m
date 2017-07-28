% OBS3_ANA  Multipath and ionospheric delay are studied as linear
%	        combinations of GPS observables. We use data from one
%	        receiver and to nine satellites.

% Written by Kai Borre, February, 23, 2002
% Copyright by Kai Borre

% Observations (in RINEX) from one satellite are read into the matrix trm1.
% The file extension is identical to the PRN. Each row corresponds
% to the observations from one epoch; the columns contain the
% following information
% 1	       2	     3	      4 	   5	    6
% P1	   Phi1      Phi2	  P2	   S1	    S2
% [m]	   [cyc]     [cyc]	  [m]	   [ ]	    [ ]
% S1 and S2 are not used; column 5 alternatively may contain P1

% Basic constants
v_light = 299792458;	 % velocity of light [m/s]
f0 = 10.23e6;		     % Basic frequency
f1 = 154*f0;		     % L1 frequency
f2 = 120*f0;		     % L2 frequency
lambda1 = v_light/f1;	 % L1 wavelength [m]
lambda2 = v_light/f2;	 % L2 wavelength [m]
alpha = (f1/f2)^2;		 % alpha

% Automatic input facility
PRN = input('Select Number of SV (1, 4, 5, 6, 8, 9, 24, 29, 30): ');
if PRN ~= [1 4 5 6 8 9 24 29 30], break, end
load(['trm1.' sprintf('%03g',PRN)]);
b = trm1;
b(:,2) = b(:,2)*lambda1; % convertion of phi [cycle] to Phi [m]
b(:,3) = b(:,3)*lambda2;

Multipath = b(:,1) + ((1+alpha)*b(:,2)-2*b(:,3))/(1-alpha);

figure(1);
plot(Multipath-mean(Multipath));
title(['Multipath for PRN ' int2str(PRN)]);
ylabel('[m]')
xlabel('Epochs [s]')

% Write code for the ionospheric delay








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end obs3_ana.m %%%%%%%%%%%%%%%%%%%
