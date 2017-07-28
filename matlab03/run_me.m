
% RUN_ME  A script to run demonstrational, basic M-code

% During my presentation in Copenhagen October 22, 2003
% I found no time to demonstrate the basic Matlab code 
% prepared for the occation.  To remedy this situation I 
% include this script which simply runs the various M-files.
% Each session is separated by a pause, so you have to press 
% enter after the end of the individual runs. Happy running!
%    
% Written October 25, 2003
%    Kai Borre

% Satconst demonstrates the GPS-satellite constallation as 
% seen in an inertial frame as well as from a rotating Earth.

satconst
close all
pause

% Via the serial i/o facility we collected data from a 
% Motorola M12 board. We stored position data as well as 
% almanac data resulting from adequate calls to the board, 
% cf. m12_data.m and m12_alm.m.

edit m12_data.tex
pause
edit m12_alm.tex
pause

% Next, the computation of a receiver position is done by 
% means of two filter versions: Kalman and Bayes. Note how the 
% PDOP value in the Kalman filter drops for each additional 
% observation filtered. This is the best illustration of the 
% filter I know of.

rec_pos
pause

% In order to make precise positioning we include phase 
% observations in addition to the basic pseudoranges. This 
% leads to an integer least-squares problem. The estimated 
% integers N are printed and the plot demonstrates the internal 
% accuracy of the observations over a time period of some 
% twenty seconds.

easy5
pause

% In the early GPS days it was common practice for planning purposes
% to compute skyplots (sterographic projections of satellite orbits)
% as well as visibility plots (when is a given satellite visible 
% with an elevation angle above, say, 15 degrees?).

skyplot
pause
close all
clear all
clc
%%%%%%%%%%%%%%%%%%%%%%% end  run_me.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
