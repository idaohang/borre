function result = dms2rad(deg,min,sec);
%DMS2RAD Conversion of degrees, minutes, and seconds to radians

%Kai Borre February 12, 1996
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date:1998/11/08  $

arg = deg+min/60+sec/3600;
result = arg*pi/180;
%%%%%%%%%%%%%%%%%%%%%%%% dms2rad.m  %%%%%%%%%%%%%%%%%%%%%%