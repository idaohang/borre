function [Phi,Lambda] = wgs2ed50(phi,lambda)
%WGS2ED50 Conversion of (phi,lambda) in WGS 84 to (Phi,Lambda)
%         in ED 50

%Kai Borre 10-17-1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date:1998/10/23  $

dtr = pi/180;
phi = input('Latitude [dd.ddddd]:  ');
lambda = input('Longitude [dd.ddddd]:  ');
phi = phi*dtr;
lambda = lambda*dtr;

a = 6378137.0;				        % Constants for WGS 84
f = 1/298.257223563;
ex = sqrt((2-f)*f)/(1-f);	     % Borre: Landmaaling (1993), (4.4)
c = a*sqrt(1+(ex)^2);
N = c/sqrt(1+(ex*cos(phi))^2);  % (4.13) and (4.16)
f1 = N*cos(phi);
x84 = f1*cos(lambda);		     % (4.46)
y84 = f1*sin(lambda);
z84 = N*(1-f)^2*sin(phi);
alpha = 0.7563093E-6;
n1 = 1+(alpha)^2;
s = 0.9999988;
x50 = s*(x84-alpha*y84)/n1+89.5;
y50 = s*(alpha*x84+y84)/n1+93.8;
z50 = s*(z84-4.5)+127.6;
Lambda = atan(y50/x50);

a = 6378388;				         % Constants for ED50
f = 1/297;
ex = sqrt((2-f)*f)/(1-f);	      % See Borre page 49
p = sqrt((x50)^2+(y50)^2);
theta = atan(z50/((1-f)*p));
t = z50+(2-f)*f*a*(sin(theta))^2*sin(theta)/(1-f);
nn = p-(2-f)*f*a*cos(theta)*(cos(theta))^2;
Phi = atan(t/nn);
N = c/sqrt(1+(ex*cos(phi))^2);
h = p/cos(phi)-N;
disp('\n Geographical coordinates transformed from WGS84 to ED50');
fprintf('\n phi_WGS84 = %10.8f and  lambda_WGS84 = %10.8f\n', ...
                              				     phi/dtr,lambda/dtr);
fprintf('\n phi_ED50  = %10.8f and  lambda_ED50  = %10.8f\n', ...
                           s						 Phi/dtr,Lambda/dtr);
%%%%%%%%%%%%%%%% end wgs2ed50.m  %%%%%%%%%%%%%%%%%%%%%%%%%













