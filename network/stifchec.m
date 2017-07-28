%STIFCHEC  Check of stiffnes matrix S for observation types 2 and 4
%          for adjacent isoscele triangles.
%          First set with SW-NE diagonal and next with NW-SE diagonal

%Kai Borre August 13, 2000
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

format rat
S1 = s2ands4(0,0,0,1,1,0,0,1); % NW SE diagonal, lower part
S2 = s2ands4(1,1,0,1,1,0,0,1); % NW SE diagonal, upper part
S3 = s2ands4(0,0,1,1,1,0,0,1); % SW NE diagonal, lower part
S4 = s2ands4(0,0,0,1,1,1,0,1); % SW NE diagonal, upper part
S1+S2;
S3+S4;
S1+S2+S3+S4
format bank
%%%%%%%%%%%%%%%%%%%%%%% end stifchec.m  %%%%%%%%%%%%%%