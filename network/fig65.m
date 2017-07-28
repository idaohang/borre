%FIG65    Difference vectors between REFDK and
%         original GI, UTM coordinates at
%         common fundamental network points

%Kai Borre February 6, 2000
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

n = 26; % # of common points
ff = fopen('diff.tex');
A = zeros(n,5);
for i = 1:n  
   line = fgetl(ff);
   A(i,1) = str2num(line(1:4));   % site
   A(i,2) = str2num(line(6:17));  % N1
   A(i,3) = str2num(line(20:32)); % E1
   A(i,4) = str2num(line(35:49)); % N2
   A(i,5) = str2num(line(52:63)); % E2
end
status = fclose(ff);
x = A(:,3);
y = A(:,2);
ze = A(:,5)-A(:,3);
zn = A(:,4)-A(:,2);
Nmax = max(y);
Nmin = min(y);
Emax = max(x);
Emin = min(x);
Ngrid = linspace(Nmin,Nmax,50);
Egrid = linspace(Emin,Emax,50);

[xi,yi] = meshgrid(Egrid,Ngrid);
zie = griddata(x,y,ze,xi,yi,'cubic'); % computing ze at gridded points
zin = griddata(x,y,zn,xi,yi,'cubic');

subplot(1,2,1); 
quiver(x,y,ze,zn)
axis([400000 700000 6050000 6450000])
axis equal, axis off

subplot(1,2,2);
hold on
[ce,he] = contour(xi,yi,zie,...
               [-.8 -.6 -.4 -.2 0.0 .2 .4 .6 .8],':r'); % Easting
[cn,hn] = contour(xi,yi,zin,...
               [1.0 1.2 1.4 1.6 1.8 2.0],'-g'); % Northing
clabel(ce,he);
clabel(cn,hn); 
axis([400000 700000 6050000 6450000])
axis equal, axis off
hold off
print -deps fig65
%%%%%%%%%%%%%%%%%%% end fig65.m  %%%%%%%%%%%%%%%%%