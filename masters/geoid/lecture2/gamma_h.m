function gamma_h = gamma_h(phi,h)
%GAMMA_H  Computes theoretical normal gravity gamma_h in m/s^2
%         for the geodetic latitude phi and at the ellipsoidal
%         height h in meters according to Somigliana's formula.
%
%Reference: Weikko A. Heiskanen & Helmut Moritz (1967): Physical
%         Geodesy. W.H.Freeman and Company

%Written by Kai Borre
%March 25, 2000

gamma_e = 9.7803253359; %m/s^2
gamma_p = 9.8321849378;
f = 1/298.257223563;
a = 6378137;
m = 0.00344978650684;

gamma = (gamma_e*(cos(phi))^2+(1-f)*gamma_p*(sin(phi))^2)/...
                   sqrt((cos(phi))^2+(1-f)^2*(sin(phi))^2);  % (2-76)
gamma_h = gamma*(1-2*(1+f+m-2*f*(sin(phi))^2)*h/a+3*(h/a)^2); % (2-123)
%%%%%%%%%%%%%%%%%%%%%%%%%% end gamma_h  %%%%%%%%%%%%%%%%%%%%%%%%
