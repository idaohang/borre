%GEO_MAP  Geoid map for the region of Denmark

%Written by Kai Borre
%March 25, 2000

%For numerical reasons it is not advisable to select
%n_max larger than 36, say. In most cases n_max =  12 is sufficient

for i = 55:57          % phi = 55:57
   for j = 9:12       % lambda = 9:12
      geoid(i,j) = wgs84_ha(i,j,0,12);
   end
end
figure;
[c,h] = contour(geoid(55:57,9:12));
clabel(c,h);
%axis([9.1 11.9 55.1 56.9])
axis equal
title('WGS84 Geoid for the Denmark Area')
%%%%%%%%%%%%%%%%%%%%%% end geo_map %%%%%%%%%%%%%%%%%%%