function fig35a(r0)
%FIG35A  Plots a Green's fuction on a unit circle 
%        with Neumann boundary conditions.
%        This Green's function is the continuous analogue
%        to the covariance matrix for 2-dimensional levelling

%Kai Borre April 24, 1995
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

G = zeros(20,20);
t = linspace(0,2*pi,20);
r = linspace(0,1,20);
for j = 1:20
  for i = 1:20
    G(i,j) = -(log(r0^2+r(i)^2-2*r0*r(i)*cos(t(j))) ...
	           +log(r(i)^2*r0^2+1-2*r(i)*r0*cos(t(j))) ...
        	       +3/2-r(i)^2-r0^2);
  end
end

h1=figure;
q = 1:1:20;
hold on
for i = 1:20
   g = G(i,:); 
   plot(q,g')
end
hold off
pause

h2 = figure;
mesh(G)
pause

h3 = figure;
figure(h3)
[X,Y,Z] = pol2cart(t,r,G);
%plot3(X,Y,Z)
surf(G)

print fig35a -deps
%%%%%%% end fig35a.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
