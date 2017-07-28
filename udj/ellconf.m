function [a,b,alpha] = ellconf(A)
% ELLCONF  Computes the eigenvalues and -vectors
%	        of a 2 x 2 covariance matrix, prints the
%	        values of the semi major, semi minor axes
%	        and the rotation angle of the pertinent
%	        confidence ellipse and plots it.

%    Written by Kai Borre
%    May 7, 1996

[m,n] = size(A);
if m ~= 2 | n ~= 2 error('Wrong dimension of matrix'); end
[V,D] = eig(A);
[lambda,k] = sort(diag(D));
V = V(:,k);
if any(any(V)) == 1, alpha = atan2(V(1,1),V(2,1)); end
rot = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
alpha = alpha*180/pi-90;  %convertion to degrees and zero to N
if alpha < 0, alpha = alpha+360; end
t = linspace(0,2*pi,49);
a = sqrt(lambda(1));
b = sqrt(lambda(2));
pl = [a*sin(t); b*cos(t)];
for t = 1:49, current = rot*pl(:,t); curve(1:2,t) = current; end

clf
hold on
axis equal
plot(curve(2,1:49),curve(1,1:49)) % interchange of x and y that
                                  % x points to N, and y to E
plot([curve(2,1),curve(2,25)],[curve(1,1),curve(1,25)],'g-')
plot([curve(2,13),curve(2,37)],[curve(1,13),curve(1,37)],'m-')
grid
hold off
print ellconf -deps
%%%%% end ellconf.m %%%%%%%%%%%%%%%
