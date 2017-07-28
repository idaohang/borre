%EASY4   The master reciever position is computed like in EASY3.
%        Next the observations taken by the rover receiver are 
%        introduced and the function baseline returns the baseline 
%        components epoch by epoch.
%	     Note that the sequence of satellites in the stored data is 
%        not the same at master and rover receivers. Therefore we 
%        must introduce a matching mechanism.

%Kai Borre 27-07-2002
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2002/07/27  $

% Read RINEX ephemerides file and convert to internal Matlab format
rinexe('SITE247J.01N','eph.dat');
Eph = get_eph('eph.dat');

% We identify the master observation file and open it
ofile1 = 'SITE247J.01O';
fid1 = fopen(ofile1,'rt');
[Obs_types1, ant_delta1, ifound_types1, eof11] = anheader(ofile1);
NoObs_types1 = size(Obs_types1,2)/2;
Pos = [];
Gdop = [];
% There are 22 epochs of data in ofile1
qend = 22;

for q = 1:qend
    [time1, dt1, sats1, eof1] = fepoch_0(fid1);
    NoSv1 = size(sats1,1);
    % We pick the observed C1 pseudoranges
    obs1 = grabdata(fid1, NoSv1, NoObs_types1);
    i = fobs_typ(Obs_types1,'C1'); 
    [pos, el, gdop] = recpo_ls(obs1(:,i),sats1,time1,Eph);
    Gdop = [Gdop gdop];
    Pos = [Pos pos];
end
me = mean(Pos,2);
spread = std(Pos,1,2)
fprintf('\nMean position as computed from %2.0f epochs:',qend)
fprintf('\n\nX: %12.3f  Y: %12.3f  Z: %12.3f\n\n', me(1,1), me(2,1), me(3,1))

figure(1);
plot((Pos(1:3,:)-Pos(1:3,1)*ones(1,q))','linewidth',2)
title(['Variation of Receiver Coordinates over ',int2str(qend),' Epochs'])
legend('X','Y','Z')
xlabel('Epochs [1 s interval]')
ylabel('[m]')
print -depsc easy4_1

% we need to close all open files and then open to read from the beginning
fclose all;

ofile1 = 'SITE247J.01O';
fid1 = fopen(ofile1,'rt');
[Obs_types1, ant_delta1, ifound_types1, eof11] = anheader(ofile1);
NoObs_types1 = size(Obs_types1,2)/2;
% Next we include the rover and identify the rover
% observation file and open it
ofile2 = 'SITE24~1.01O';
fid2 = fopen(ofile2,'rt');
[Obs_types2, ant_delta2, ifound_types2, eof12] = anheader(ofile2);
NoObs_types2 = size(Obs_types2,2)/2;
master = me;  % best possible estimate of master position
bases = [];

for q = 1:qend
    [time1, dt1, sats1, eof1] = fepoch_0(fid1);
    [time2, dt2, sats2, eof2] = fepoch_0(fid2);
    if time1 ~= time2
        disp('Epochs not corresponding')
        break
    end;
    NoSv1 = size(sats1,1);
    NoSv2 = size(sats2,1);
    % We pick the observations
    obsm = grabdata(fid1, NoSv1, NoObs_types1);
    obsr = grabdata(fid2, NoSv2, NoObs_types2);
    i = fobs_typ(Obs_types1,'C1');
    obs1 = obsm(:,i);
    obs2 = obsr(:,1);
    for s = 1:NoSv1
        ind = find(sats1(s) == sats2(:));
        obs(s,1) = obs1(s)-obs2(ind);
    end
    base = baseline(master,obs,sats1,time1,Eph);
    bases = [bases base];
end
me1 = mean(bases,2);
spread1 = std(bases,1,2)
fprintf('\nBaseline Components as Computed From %2.0f Epochs:',qend)
fprintf('\n\nX: %12.3f  Y: %12.3f  Z: %12.3f', me1(1,1),me1(2,1),me1(3,1))

figure(2);
plot((bases(1:3,:)-bases(1:3,1)*ones(1,q))','linewidth',2)
title(['Variation of Baseline Coordinates over ',int2str(qend),' Epochs'])
legend('X','Y','Z')
xlabel('Epochs [1 s interval]')
ylabel('[m]')
print -depsc easy4_2

figure(3);
plot(Gdop,'linewidth',2)
axis([1 length(Gdop) 0 5])
title('GDOP')

%%%%%%%%%%%%%%%%%%% end easy4.m %%%%%%%%%%%%%%%



