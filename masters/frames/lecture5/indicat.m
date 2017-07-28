%INDICAT   Plot of indicatrices for part of the UTM mapping

%Written by Kai Borre
%April, 21, 2000

dtr = 180/pi; % degree to radian
figure;
hold on
for phi = 0:20:80
   for lambda = -6:3:6
      [N,E] = geo2utm(phi,lambda,30); % call in degrees
      gamma = sin(phi/dtr)*lambda/dtr;
      m = 1+E^2/(2*6370000^2);
      scale = 1.e5; % scale to make the indicatrix visible
      draw_ell(E,N,[(scale*m)^2 0;0 (scale*m)^2]);
      deltaE = scale*gamma;
      deltaN = scale;
      plot([E E-3*deltaE],[N N+3*deltaN])
   end
end
title('Indicatrices for the UTM mapping')
hold off
%%%%%%%%%%%%%%%%%%%%%% end indicat.m  %%%%%%%%%%%%%%%