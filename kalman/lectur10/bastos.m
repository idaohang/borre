%BASTOS  Test of GPS code and phase observations as described in
%           BASTOS L. Bastos and H. Landau (1988)
%	           Fixing cycle slips in dual-frequency
%	           kinematic GPS-applications using Kalman filtering.
%	           manuscripta geodetica, volume 13, pages 249--256

%Kai Borre, February 3, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date:1999/01/02 $

% Automatic input facility
receiver = input('Select master or rover (m, r): ','s');
if receiver ~= ['m' 'r'], break, end
sv = input('Select PRN (2, 9, 16, 23, 26, 27): ');
if sv ~= [2 9 16 23 26 27], break, end
datasv = ['one_' receiver int2str(sv)];
filename = [datasv '.dat'];
fid = fopen(filename);
data = fread(fid,inf,'double');
r = size(data,1);
B = reshape(data,r/5,5);
range = B(5:r/5,1:4);
elevation = B(5:r/5,5);

%Repair of clock reset; only for pseudoranges
spikes = diff(range(:,1));
i = find(abs(spikes) > 280000);
for j = 1:size(i)
   if spikes(i) < 0
      corr = 299792.458; 
   else
      corr = -299792.458;
   end
   for k = i(j)+1:size(range,1)
      range(k,1) = range(k,1)+corr;
      range(k,3) = range(k,3)+corr;
   end
end

r = size(range,1);
rt = 1:r;

% P1(m)   Phi1(m)    P2(m)	   Phi2(m)    elevation angle(degree)
delta = (range(:,2)-range(1,2))-(range(:,1)-range(1,1)); % Eq. (2-1)
% Ionospheric delay as derived from code on L1 and L2
iono = (range(:,3)-range(1,3))-(range(:,1)-range(1,1));
% If no cycle slips occur delta_phi is only influenced by
% ionospheric effects; range(:,4) is already multiplied by f1/f2!
delta_phi = range(:,2)-range(1,2)-(range(:,4)-range(1,4)); % Eq. (3-1)

%Initialization of Kalman filter
Delta_t = 1;
%F = eye(3)+[0 1 0;0 0 1;0 0 0]*Delta_t+ ...
%		                  [0 0 1;0 0 0;0 0 0]*Delta_t^2/2+...
F =[1 Delta_t Delta_t^2/2;0 1 Delta_t;0 0 1];
Q = diag([1.e-6 1.e-8 1.e-10]);
P = 0.1*eye(3);
R = 0.0025;
x = zeros(3,1);
A = [1 0 0];
delta_filt = [];

for i = 1:r
   [x,P,c1,c2] = k_updatf(x,P,A,delta(i),R,Q,F);
   delta_filt = [delta_filt x];
end

fig1 = figure;
set(fig1,'Numbertitle','off','name',...
   ['Different combinations of GPS-observations for PRN ' int2str(sv)]);
ax1 = gca;
line(flipud(rt+3),flipud(elevation),'color','r','parent',ax1);
set(ax1,'YAxisLocation','right',...
   'Color','none',...
   'xtick',[],'XTickLabel',[],...
   'XLim',get(ax1,'XLim'),...
   'YLim',[0 90],...
   'fontsize',8,...
   'xdir','reverse')
set(get(gca,'ylabel'),'String',('Elevation Angle'),'fontsize',8,...
   'color','r','verticalalignment','bottom','rotation',270)
ax2 = axes('pos',get(ax1,'pos'),'fontsize',8);
hold on
line(rt,delta,'color','b','linestyle','--','parent',ax2);
line(rt,delta_filt(1,:),'color','k','linestyle','-.','parent',ax2);
line(rt,iono,'color','g','linestyle',':','parent',ax2);
line(rt,delta_phi,'LineWidth',2,'parent',ax2);
legend('Code-Phase combination','Filtered Code-Phase combination',...
           'Ionospheric delay','Phase combination');
legend('Code-Phase combination','Filtered Code-Phase comb.',...
           'Ionospheric delay','Phase combination');
xlabel('Epochs (20 s interval)')
ylabel('Meters')
title(...
   ['Different combinations of GPS-observations for PRN ' int2str(sv)]);
axes(ax1); % to make elevation line visible
hold off

print -deps bastos
%%%%%%%%%%%% end bastos.m  %%%%%%%%%%%%%%%%%%%%
