%RAWRECPO  Plot of raw 3-D receiver positions as created 
%          by a call like bdata_po('b0005a94.076').
%          The data are stored in bdata_*.dat

%Written by Kai Borre, 30 April 1999; 
%revised August 21, 1999

fidm = fopen('bdata_ma.dat'); 
[pos, nums] = fread(fidm,Inf,'double');
Posm = reshape(pos,nums/3,3);
Posm(1:60,:) = []; % master receiver starts 30 seconds ahead of rover
MPosm = mean(Posm);
[pm qm] = size(Posm);

fidr = fopen('bdata_ro.dat'); 
[pos, nums] = fread(fidr,Inf,'double');
Posr = reshape(pos,nums/3,3);
MPosr = mean(Posr);
[pr qr] = size(Posr);

epochs = min(pm,pr);
M = Posm(1:epochs,:)-ones(epochs,1)*MPosm;
R = Posr(1:epochs,:)-ones(epochs,1)*MPosr;

figure(1);
subplot(2,1,1), plot(1:epochs,M)
set(get(gca,'Title'),'String',['\it{X, Y, Z}\rm components',...
            ' of receiver position relative to mean position'])
ylabel('[m]')
xlabel('[s]')

for i = 1:epochs
   [dphi(i),dlambda(i)] = togeod(6378137,298.257223563, ...
                             Posm(i,1),Posm(i,2),Posm(i,3));
end;

subplot(2,1,2), plot(dlambda,dphi)
title('Raw receiver positions relative to mean position')
set(get(gca,'YLabel'),'String','\it{\phi}')
set(get(gca,'XLabel'),'String','\it{\lambda}')
axis equal

figure(2);
start = 300; % the first 300 epochs are not good
dpos = R(start:epochs,:)-M(start:epochs,:);
subplot(2,1,1), plot(start:epochs,dpos) 
set(get(gca,'Title'),'String',['\it{X, Y, Z}\rm  components', ...
                                       ' of differential vector'])
ylabel('[m]')
xlabel('[m]')

[phi,lambda] = togeod(6378137,298.257223563,...
                            MPosm(1),MPosm(2),MPosm(3));
for i = 1:epochs-start
   [Az(i),El(i),dist(i)] = topocent(MPosm,dpos(i,:)');
end
[x,y] = pol2cart(Az*pi/180,dist);
subplot(2,1,2), plot(y-y(1),x-x(1))
title('Differential receiver position')
ylabel('[m]')
xlabel('[m]')
axis equal

%%%%%%%%%%%%%%%%% rawrecpo.m  %%%%%%%%%%%%%%%%%%