%COPY12  The master reciever position is computed like in EASY3.
%        Next the observations taken by the rover receiver are 
%        introduced and the function baseline returns the baseline 
%        components epoch by epoch.
%	     Note that the sequence of satellites in the stored data is 
%        not the same at master and rover receivers. Therefore we 
%        must introduce a matching mechanism.

%Kai Borre 25-11-2002
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2002/07/27  $

% Read RINEX ephemerides file and convert to internal Matlab format
rinexe('SITE247J.01N','eph.dat');
Eph = get_eph('eph.dat');

% We identify the master observation file and open it
ofile1 = 'SITE247j.01O';
fid1 = fopen(ofile1,'rt');
[Obs_types1, ant_delta1, ifound_types1, eof11] = anheader(ofile1);
NoObs_types1 = size(Obs_types1,2)/2;
% There are 22 epochs of data in ofile1
qend = 22;

for q = 1:qend
    [time1, dt1, sats1, eof1] = fepoch_0(fid1);
    NoSv1 = size(sats1,1);
    % We pick the observed C1 pseudoranges
    obs1 = grabdata(fid1, NoSv1, NoObs_types1);
    i = fobs_typ(Obs_types1,'C1'); 
    [pos(:,q), el, gdop] = recpo_ls(obs1(:,i),sats1,time1,Eph);
end
me = mean(pos,2);
[phi_i,lambda_i,h_i] = togeod(6378137,298.257223563,me(1),me(2),me(3));
spread = std(pos,1,2)
fprintf('\nMean position as computed from %2.0f epochs:',qend)
for i = 1:qend
    [e(i),n(i),u(i)] = xyz2enu(phi_i,lambda_i,pos(1,i)-me(1),pos(2,i)-me(2),pos(3,i)-me(3));
end
fprintf('\n\nX: %12.3f  Y: %12.3f  Z: %12.3f\n\n', me(1,1), me(2,1), me(3,1))
figure(1);
plot(1:qend,[(e-e(1))' (n-n(1))' (u-u(1))'],'linewidth',2)
title(['Variation of Receiver Coordinates over ',int2str(qend),' Epochs'])
legend('E','N','U')
xlabel('Epochs [1 s interval]')
ylabel('[m]')
print -depsc easy4_1
fclose all;
% We identify the rover observation file and open it
ofile1 = 'SITE24~1.01O';
fid1 = fopen(ofile1,'rt');
[Obs_types1, ant_delta1, ifound_types1, eof11] = anheader(ofile1);
NoObs_types1 = size(Obs_types1,2)/2;
Posr = [];
% There are 22 epochs of data in ofile2
for q = 1:qend
    [time1, dt1, sats1, eof1] = fepoch_0(fid1);
    NoSv1 = size(sats1,1);
    % We pick the observed C1 pseudoranges
    obs1 = grabdata(fid1, NoSv1, NoObs_types1);
    i = fobs_typ(Obs_types1,'C1'); 
    [posr(:,q), el, gdop] = recpo_ls(obs1(:,i),sats1,time1,Eph);
end
mer = mean(posr,2);
spread = std(posr,1,2)
fprintf('\nMean position as computed from %2.0f epochs:',qend)
for i = 1:qend
    [er(i),nr(i),ur(i)] = xyz2enu(phi_i,lambda_i,posr(1,i)-mer(1),posr(2,i)-mer(2), ....
        posr(3,i)-mer(3));
end

fprintf('\n\nX: %12.3f  Y: %12.3f  Z: %12.3f\n\n', mer(1,1), mer(2,1), mer(3,1))

figure(2);
plot(1:qend,[(er-er(1))' (nr-nr(1))' (ur-ur(1))'],'linewidth',2)
title(['Variation of Receiver Coordinates over ',int2str(qend),' Epochs'])
legend('E','N','U')
xlabel('Epochs [1 s interval]')
ylabel('[m]')

figure(3);
plot(1:qend,[(e-er)' (n-nr)' (u-ur)'],'linewidth',2)
break
plot3(Pos(1,:)-Posr(1,:),Pos(2,:)-Posr(2,:),Pos(3,:)-Posr(3,:),'o')
me-mer
print -depsc easy4_1
%%%%%%%%%%%%%%%%%%% end easy4.m %%%%%%%%%%%%%%%

