% EASY8    Test for cycle slip and repair of receiver clock offset. 
%          We use the ionosphere free linear combination.
%          Formulas are described in Strang and Borre, page 491

%Kai Borre 27-07-2002
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2002/07/27  $

v_light = 299792458;    % vacuum speed of light m/s
f0 = 10.23*10^6;
f1 = 154*f0;
f2 = 120*f0;
lambda1 = v_light/f1;
lambda2 = v_light/f2;
alpha1 = f1^2/(f1^2-f2^2);
alpha2 = 1-alpha1;

ofile = 'site247j.01o';
%ofile = 'kofi2.01o';

fid = fopen(ofile,'rt');
[Obs_types, ant_delta, ifound_types, eof1] = anheader(ofile);
if ((ifound_types == 0) | (eof1 == 1))
    error('Basic information is missing in RINEX file')
end;
NoObs_types = size(Obs_types,2)/2;
j = fobs_typ(Obs_types,'P1');
k = fobs_typ(Obs_types,'P2');
l = fobs_typ(Obs_types,'L1');
m = fobs_typ(Obs_types,'L2');
cols = [j k l m];
fid = fopen(ofile,'rt');

for i = 1:20
    [tr_RAW, dt, sv, eof2] = fepoch_0(fid); 
    NoSv = size(sv,1);
    obs = grabdata(fid, NoSv, NoObs_types);  
    Obs = obs(:,cols);
    for j = 1:NoSv
        if i == 1
            Obs0 = Obs;
            P(:,1) = alpha1*Obs0(:,1)+alpha2*Obs0(:,2);
            Phi(:,1) = alpha1*lambda1*Obs0(:,3)+alpha2*lambda2*Obs0(:,4);
        else
            P(j,i) = alpha1*Obs(j,1)+alpha2*Obs(j,2);
            deltaP(j,i) = P(j,i)-P(j,1);
            Phi(j,i) = alpha1*lambda1*Obs(j,3)+alpha2*lambda2*Obs(j,4);
            deltaPhi(j,i) = Phi(j,i)-Phi(j,1);
        end
    end
end
fclose all;
DP = diff(deltaP,1,2);
% Repair of clock reset of 1ms ~ 299 km; affects only pseudoranges
i1 = find(abs(DP(1,:)) > 280000);

for j = i1
   if DP(:,j) < 0
       DP(:,j) = DP(:,j)+299792.458; 
   else 
       DP(:,j) = DP(:,j)-299792.458; 
   end
end
DPhi = diff(deltaPhi,1,2);
DDP = diff(DP,1,2);
DDPhi = diff(DPhi,1,2);

figure(1);
subplot(2,2,1), plot(DP')
set(get(gca,'ylabel'),'String',({'Pseudoranges Differenced','over Time [m]'}),...
    'fontsize',12,'color','r','verticalalignment','bottom','rotation',90)   
subplot(2,2,2), plot(DPhi')
set(get(gca,'ylabel'),'String',({'Phases Differenced','over Time [m]'}),...
    'fontsize',12,'color','r','verticalalignment','bottom','rotation',90)   
subplot(2,2,3), plot(DDP')
set(get(gca,'ylabel'),'String',({'Double Differenced','Pseudoranges [m]'}),...
    'fontsize',12,'color','r','verticalalignment','bottom','rotation',90)   
subplot(2,2,4), plot(DDPhi')
set(get(gca,'ylabel'),'String',({'Double Differenced','Phases [m]'}),...
    'fontsize',12,'color','r','verticalalignment','bottom','rotation',90)   

%ax = axes('units','normal','position',[.075 .075 .9 .9],'visible','off');
%set(get(ax,'title'),'visible','on')
%title('Check of Cycle Slips','fontsize',16,'fontweight','bold')
%h = get(ax,'title');

print -deps easy8
%%%%%%%%% end easy8.m %%%%%%%%%
