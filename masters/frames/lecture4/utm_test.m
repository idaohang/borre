% Example from UTM manual
% hypothetical station in Algeria
% International ellipsoid

% Written by Kai Bore
% March 31, 2000

phi = dms2rad(34,15,34.742)*180/pi;
lambda = dms2rad(5,57,16.842)*180/pi;
[N,E] = geo2utm(phi,lambda,31);

% Reverse computation

[b,l] = utm2geo(N,E,31);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end utm_test.m %%%