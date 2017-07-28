%Script for Lecture 2

%Written by Kai Borre
%October 29,1998

%We want to compute (X,Y,Z) for all available satellites
%at a given epoch. 

%We start from an e-file created by hosing a Z12-receiver.
%It was named e0810a94.076. This means that the file is created
%in year 1994 on day 76 = March 17.
%We want to changed the format of this binary file into a 
%matrix of ephemerides. Each ephemeris is contained in a column.
%All columns have 21 rows. The first row contains the PRN
%number, etc. For a full description of the ephemeris,
%see the edata.m file.

edata('e0810a94.076','edata.dat');

%We read the ephemerides to Eph

Eph = get_eph('edata.dat');

%Next we have to specify the PRN's for which we want to compute
%the position and the time. We get the number of all 
%available satellites by the following call

available = Eph(1,:)

%Immediately you see that the individual PRN's are repeated a lot
%of times. We only want a list of the PRN's

PRNS = unique(Eph(1,:))
NO_PRN = size(PRNS,2)

%Now we know all PRNs above the local horizon. Next we must select
%the time for computing the position of the satellites.
%That time is prescribed by the epoch when the observation was made.
%These epochs are stored in the b-file.

%In the following we work specifically with file 00050761.94o.
%This file may be created by means of the ashtorin file in GPPS.
%However we have used another formatter. We focus on the last 
%epoch in the file. It is tagged 
%94  3 17 10 12 0.00
%We need to convert this epoch time into seconds of week.

jd = julday(1994,3,17,10.2);
[week,SOW] = gps_time(jd)

%Finally we compute the satellite positions at time SOW

for i = 1:NO_PRN
   col_Eph = find_eph(Eph,PRNS(i),SOW);
   X(i,:) = satpos(SOW,Eph(:,col_Eph))';
end

%Now we have computed the positions of all available PRN's.
%We store the PRN's and the coordinates X in 

B = [PRNS' X]

%%%%%%%%%%%%%%%%%%%%%%%% end run.m  %%%%%%%%%%%%%%%%%



