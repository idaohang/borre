%SAV	Plots of a satellite's acceleration and velocity

%Kai Borre, October 28, 1997
%Copyright (c) by Kai Borre
%$Revision: 1.0 $ $Date:1998/11/09  $

eph = get_eph('pta.nav');
for t = 1:100
   x(t,:) = satpos(170000+t,eph(:,1))';
end
figure
plot(diff(diff(x)));
xlabel('Seconds');
ylabel('X, Y, and Z components of Acceleration [m/s^2]');

for t = 1:99
   dist(t) = norm(x(t,:)-x(t+1,:));
end
figure
plot(dist)
xlabel('Seconds');
ylabel('Velocity of satellite [m/s]');
%%%%%%%%% end sav.m  %%%%%%%%%%%%%%%%%
