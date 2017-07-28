%DK97GEOI  Interpolation of geoid undulation in Denmark
%          Plot of geoid contours

% Written by Kai Borre
% April 23, 2000

%%fid = fopen('refdk97.n','r');
%%G = fscanf(fid,'%g %g %g %g %g',[5 inf]); 
load G
%G = G';
[m,n] = size(G);
%Phi = G(:,2);
%Lambda = G(:,3);
%N = G(:,5);
figure(1);
axis([8 16 54 58])
hold on
plot(G(:,3),G(:,2),'+')
title(['DKGEOID97 and ' int2str(m) ' REFDK-Points (+)'])
xlabel('Longitude [\circ]')
ylabel('Latitude [\circ]')
ml = 200;
nl = 200;
ylin = linspace(54,58,nl); % we interpolate between phi = 54 and
                           % 58, longitude lambda = 8 and 16
xlin = linspace(8,16,ml);
[X,Y] = meshgrid(xlin,ylin);
Z = griddata(G(:,3),G(:,2),G(:,5),X,Y,'cubic');
v = [34:.1:41];
[C,h] = contour(X,Y,Z,v);
clabel(C,h,[33 34 35 36 37 38 39 40 41]); 
grid on
hold off
%%%%%%%%%%%%%%%%%%%%% end dk97geoi.m  %%%%%%%%%%%%%%%%%%