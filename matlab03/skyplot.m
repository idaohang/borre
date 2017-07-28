function skyplot(almanac,mask,phi,lambda)
%          Makes a stereographic plot of GPS satellite orbits from an almanac.
%          The plot is as seen from the position (phi,lambda). All orbits with
%          elevation angle lower than mask are omitted.

% Written by Kai Borre
% November 8, 2002

% typical call
% skyplot('alm.dat',10,[67 0 0],[-50 0 0])

if nargin == 0
    almanac = 'alm.dat';
    mask = 10;
    phi = [57 0 0];
    lambda = [10 0 0];
end

%reading the ephemerides
fide = fopen(almanac,'r');
Eph = fread(fide,inf,'double');
m = length(Eph);
eph = reshape(Eph,21,m/21);

%transformation of location coordinates (phi,lambda,h) to (X,Y,Z)
Phi = dms2rad(phi(1),phi(2),phi(3));
Phi = Phi*180/pi;
Lambda = dms2rad(lambda(1),lambda(2),lambda(3));
Lambda = Lambda*180/pi;
[M(1,1),M(2,1),M(3,1)] = frgeod(6378137,298.257222101,Phi,Lambda,0);

figure(1);
set(gcf,'UserData',zeros(2,1)*inf);
plot1 = polarhg(0,0,'.'); 

satsum = zeros(1,288);
visible = zeros(2*size(eph,2)+1,288);

for i = 1:size(eph,2)
    Az = [];
    El = [];
    for j = 1:288 % step size for j: 300 s = 5 min
        time = 300*(j-1);
        S = satpos(time,eph(:,i));
        [az,el,d] = topocent(M,S-M);
        if el > mask
            Az = [Az az];
            El = [El el];
            satsum(j) = satsum(j)+1;
            visible(2*i,j) = 1;
        end
    end
    polarhg(Az*pi/180,90-El,'tdir','clockwise','rlim',[0 90],'rtick',[0 30 60 90],... 
        'torig','up','color','b','linestyle','.')
    pause(.5)  % rtick 90 60 30 0
    hold on
end    
title({['Skyplot for the position (\phi,\lambda) = (' num2str(phi(1)) '\circ,' ...
            num2str(lambda(1)) '\circ)'],['Elevation mask ' num2str(mask) '\circ' ]})
text(-120,-160,['PRNs: ' num2str(eph(1,:))])
print -depsc skyplot1

figure(2);
fig21 = subplot(2,1,1);
area(satsum)
set(fig21,'XTick',1:71:288)
set(fig21,'XTickLabel',{'0','6','12','18','24'})
ylabel('Number of Visible Satellites')
title(['Elevation mask ' num2str(mask) '\circ'])

fig22 = subplot(2,1,2);
imagesc(flipud(visible)); %colormap(gray)
set(fig22,'XTick',1:71:288)
set(fig22,'XTickLabel',{'0','6','12','18','24'})
set(fig22,'YTickLabel',[])
title({'The Horizontal Lines Represent Satellites,','Increasing Numbers from Top to Bottom'})
xlabel('Red Color Indicates Satellite Up')
print -depsc skyplot2

%%%%%%%%%%%%%%%%%%%%% end skyplot.m %%%%%%%%%%%%%%
