% NORMAT  Script to compute GPS receiver position using various
%         methods
%         1) iterative least-sqaures method
%         2) a least-squares search technique
%         3) Bancroft's algorithm
%         4) Kleusberg's algorithm
%         5) Bayes filter
%         6) Kalman filter

%Written by Kai Borre July 13, 2002
%Copyright by Kai Borre

% Our computational example is based upon data collected by a 
% JPS Eurocard receiver on April 25, 2001 at Aalborg, Denmark.

v_light = 299792458;
% The epoch time in seconds of week, sow: 
time = 207600;                    
% The receiver tracked the following PRNs
prns = [1; 4; 7; 13; 20; 24; 25]; 
% The observed pseudoranges at the epoch time were
P = ...
    [20532016.23;  
    21255529.02;
    24648801.64;
    21267722.01;
    21900015.95;
    23828511.82;
    24104654.81];
% The ephemerides are needed for computation of the satellite position at 
% the transmit time of the signal. That time is the receiver epoch time minus 
% the travel time \tau^k. \tau^k is the measured pseudorange P divided  
% by the velocity of light.

Eph = get_eph('eph.dat');

% 1) the variable B contains the ECEF coordinates of the PRN's as 
%    well as pseudoranges corrected for satellite clock offsets and tropospheric 
%    delay. El is the elevation angle of the  
%    PRNs and GDOP tells about the general quality of the actual satellite 
%    configuration
[pos_ls, El, GDOP, B] = recpo_ls(P,prns,time,Eph); 
fprintf('\nOrdinary Least-Squares Solution\n')
fprintf('\n X= %6.2f [m]  Y= %6.2f [m] Z= %6.2f [m] c dt= %6.2f [m]\n', ...
                pos_ls(1),pos_ls(2),pos_ls(3),pos_ls(4))

% 2) We start making a call to togeod to get a value for the height h
[p,l,h] = togeod(6378137,298.257223563,pos_ls(1),pos_ls(2),pos_ls(3));
[phi,lambda] = recpos_s(B);
[X,Y,Z] = frgeod(6378137,298.257223563,phi,lambda,h);
fprintf('\nResult of Search Technique\n')
fprintf('\n X= %6.2f [m]  Y= %6.2f [m] Z= %6.2f [m]\n', X,Y,Z);

% 3) 
pos_ban = bancroft(B);
fprintf('\nBancroft Algorithm\n')
fprintf('\n X= %6.2f [m]  Y= %6.2f [m] Z= %6.2f [m] c dt= %6.2f [m]\n', ...
                pos_ban(1),pos_ban(2),pos_ban(3),pos_ban(4));

% 4)
pos_kle = kleus(B(1:4,:));
fprintf('\nKleusberg Algorithm\n')
fprintf('\n X= %6.2f [m]  Y= %6.2f [m] Z= %6.2f [m]\n', ...
                pos_kle(1),pos_kle(2),pos_kle(3));

% 5)
pos_bayes = b_point(Eph,P,prns,time); 
fprintf('\nBayes Filter\n')
fprintf('\n X= %6.2f [m]  Y= %6.2f [m] Z= %6.2f [m] c dt= %6.2f [m]\n', ...
                pos_bayes(1),pos_bayes(2),pos_bayes(3), pos_bayes(4)*v_light*1.e-9);

% 6)
pos_kalman = k_point(Eph,pos_bayes,P,prns,time);
fprintf('\nKalman Filter\n')
fprintf('\n X= %6.2f [m]  Y= %6.2f [m] Z= %6.2f [m] c dt= %6.2f [m]\n', ...
                pos_kalman(1),pos_kalman(2),pos_kalman(3),pos_kalman(4)*v_light*1.e-9);

%%%%%%%%%%%%%%%%%%%%%%%% end normat.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%