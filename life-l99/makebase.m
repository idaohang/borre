%MAKEBASE  We assume that you already have run the M-file bdata 
%          with b-files from master and rover. The call would be
%             bdata('masterfile','roverfile') or specifically
%	           bdata('b0810a94.076','b0005a94.076')
%          That creates the file bdata.dat.
%          Likewise we assume that the calls
%             edata('masterfile','roverfile') and
%             sdata('masterfile','roverfile')
%          have created edata.dat and sdata.dat.
%          Making Double Differences of Code and Phase Observations.
%          Estimates ambiguities on L1 and L2, and baseline components.
%
%          References refer to
%          Strang, Gilbert and Borre, Kai (1997): Linear Algebra,
%          Geodesy, and GPS. Wellesley-Cambridge Press.

%Kai Borre, August 25, 1999

%THIS SAMPLE CODE DOES NOT ACCOUNT FOR CYCLE SLIPS

% Initial computations of constants
v_light = 299792458; 	     % vacuum speed of light m/s
f1 = 154*10.23E6;		        % L1 frequency Hz
f2 = 120*10.23E6;			     % L2 frequency Hz
lambda1 = v_light/f1;	     % wavelength on L1:  .19029367  m
lambda2 = v_light/f2;	     % wavelength on L2:  .244210213 m

tic

fidb = fopen('bdata.dat');
[da,count] = fread(fidb,Inf,'double');
rows = count/7;
B = reshape(da,rows,7);
clear da
i1 = [];
for i = 1:rows-1
   if B(i+1,1) < B(i,1)
      i1 = [i1 i];
   end;
end
BR = B(1:i1,:);	         % rover data
BM = B(i1+1:rows,:);       % master data
clear B

% Further investigation of master data
prns = unique(BM(:,2));
[mr,mc] = size(BM);
time_beg = BM(1,1);
time_end = BM(mr,1);
fprintf('\nData for master')
fprintf('\nBegin time %10.0f sow', time_beg)
fprintf('\nEnd time   %10.0f sow', time_end)
noepo = find(BM(:,1) > BM(1,1));
epoch_interval = BM(noepo(1),1)-BM(1,1);
fprintf('\nEpoch interval   %4.0f sec\n', epoch_interval)
epochs = (time_end-time_beg)/epoch_interval+1;
%We multiply by Inf to get a matrix of NaN's
%Each row in an obs matrix has the following five columns   
%  P1	  Phi1   P2    Phi2   elevation
obs = zeros(epochs,5*length(prns))*Inf; 

for i = 1:length(prns)       %runs over prns
   for j = 1:mr              %runs over all rows in BM
      if prns(i,1) == BM(j,2)
         k = (BM(j,1)-BM(1,1))/epoch_interval + 1;
         obs(k,(i-1)*5+1:(i-1)*5+5) = BM(j,3:7);            
      end
   end   
end

%Further investigation of rover data
prnsr = unique(BR(:,2));
[rr,rc] = size(BR);
time_beg = BR(1,1);
time_end = BR(rr,1);
fprintf('\nData for rover')
fprintf('\nBegin time %10.0f sow', time_beg)
fprintf('\nEnd time   %10.0f sow', time_end)
noepor = find(BR(:,1) > BR(1,1));
epoch_interval = BR(noepor(1),1)-BR(1,1);
fprintf('\nEpoch interval   %4.0f sec\n', epoch_interval)
epochs = (time_end-time_beg)/epoch_interval+1;
obsr = zeros(epochs,5*length(prns))*Inf;

for i = 1:length(prns)       %runs over prns
   for j = 1:rr              %runs over all rows in BR
      if prns(i,1) == BR(j,2)
         k = (BR(j,1)-BR(1,1))/epoch_interval + 1;
         obsr(k,(i-1)*5+1:(i-1)*5+5) = BR(j,3:7);            
      end
   end   
end

%Finding the PRN with largest mean elevation at master
ME = mean(obs(:,5*(1:length(prns)))); 
[y,refprn] = max(ME);
%We delete all satellites with mean elevation lower than 15 degrees
kprn = [1 2 4 6 7]; %%find(ME > 15);
ttprn = kprn;
kprn = setdiff(kprn,refprn);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start of ambiguity estimation, see pages 489--490
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n] = size(obs(:,1));
% Definition of filter matrix, see equation (15.9)
Ad = [1 0       0;
      1 lambda1 0;
      1 0       0;
      1 0       lambda2];
sd = [0.3 0.005 0.3 0.005]; % standard deviation of observations
W = diag(1./sd.^2);         % diagonal weight matrix
AW = Ad'*W;
N = AW*Ad;
ef = 219;   % first epoch
el = 718;   % last epoch
rk = 0;
x_amb = [];

for k = kprn
   %The observation vector b contains double differences 
   %of     P1   Phi1   P2    Phi2. All data in a row pertain to
   %the same epoch
   b = [obs(:,(refprn-1)*5+1)-obs(:,(k-1)*5+1)-...
         obsr(:,(refprn-1)*5+1)+obsr(:,(k-1)*5+1), ...
         (obs(:,(refprn-1)*5+2)-obs(:,(k-1)*5+2)-...
         obsr(:,(refprn-1)*5+2)+obsr(:,(k-1)*5+2))*lambda1, ...
         obs(:,(refprn-1)*5+3)-obs(:,(k-1)*5+3)-...
         obsr(:,(refprn-1)*5+3)+obsr(:,(k-1)*5+3), ...
         (obs(:,(refprn-1)*5+4)-obs(:,(k-1)*5+4)-...
         obsr(:,(refprn-1)*5+4)+obsr(:,(k-1)*5+4))*lambda2]; 
   N22 = zeros(2,2);
   RS21 = zeros(2,1);
   % i runs over epochs
   for i = 300:400 %selected to get reliable ambiguities  
      RS = AW*b(i,:)';
      N22 = N22+N(2:3,2:3)-N(2:3,1)*N(2:3,1)'/N(1,1);
      RS21 = RS21+RS(2:3,1)-N(2:3,1)*RS(1,1)/N(1,1);
   end %i
   x = inv(N22)*RS21;
   x_amb = [x_amb x];
   K1 = round(x(1)-x(2));
   K2 = round(60*x(1)-77*x(2));
   trueN2 = round((60*K1-K2)/17);
   trueN1 = round(trueN2+K1);
   fprintf(['\nAmbiguities for PRNs'...
         '%3.0f and %2.0f'],prns(refprn),prns(k));     
   fprintf('\nResult: N_L1: %10.0f, N_Lw: %10.0f',...
      trueN1, trueN1-trueN2); 
   rk = rk+1;
   amb(rk,:) = [trueN1 trueN2];
end
fprintf('\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATION OF VECTOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%master position
pos = [3427907.43;
        603542.18;
       5326902.62];
X_i = pos(1:3);
x = zeros(3,1); %first guess for preliminary vector components
X_j = X_i+x;
[phi_i,lambda_i,h_i] = ...
   togeod(6378137,298.257223563,X_i(1),X_i(2),X_i(3));

% Calculation of weight matrix
noprn = length(kprn);
D = [ones(noprn,1) -eye(noprn) -ones(noprn,1) eye(noprn)];
C = inv(D*D');

% We relate each satellite with a fixed ephemeris
Eph = get_eph('edata.dat');
time = BR(1,1);
for t = 1:noprn
   col_Eph(t) = find_eph(Eph,prns(kprn(t)),time);
end
% and the reference satellite
col_Eph_r = find_eph(Eph,prns(refprn),time);

for slides = 1:3
   switch slides
   case 1
      Amb = x_amb(1,:)';
   case 2, 3
         Amb = amb(:,1); 
   end;   
   ressum = zeros(noprn,1)*inf;
   Var = 100;  % Initial value for variance of an epoch solution
   
   [phi_j,lambda_j,h_j] =  ...
      togeod(6378137,298.257223563,X_j(1),X_j(2),X_j(3));
   Normal = zeros(3,3);
   RightSide = zeros(3,1);
   xx = [];
   varsum = 0;
   t1 = 0;
   
   for p = ef:el
      if slides == 3 & p == 469, Amb(1) = Amb(1)-1; end
      time = BM(1,1)+(p-1)*epoch_interval;
      [rhok_j,Xk_ECF] = get_rho(time, ...
         obsr(p,(refprn-1)*5+1), Eph(:,col_Eph_r), X_j);
      [rhok_i,Xk_ECF] = get_rho(time, ...
         obs(p,(refprn-1)*5+1), Eph(:,col_Eph_r), X_i);
      tt = 0;
      for t = kprn
         tt = tt+1;
         [rhol_j,Xl_ECF] = get_rho(time,...
            obsr(p,(t-1)*5+1), Eph(:,col_Eph(tt)), X_j);
         [rhol_i,Xl_ECF] = get_rho(time, ...
            obs(p,(t-1)*5+1), Eph(:,col_Eph(tt)), X_i);
         A(tt,:) = [(Xk_ECF(1)-X_j(1))/rhok_j ...
               - (Xl_ECF(1)-X_j(1))/rhol_j  ...
               (Xk_ECF(2)-X_j(2))/rhok_j - (Xl_ECF(2)-X_j(2))/rhol_j ...
               (Xk_ECF(3)-X_j(3))/rhok_j - (Xl_ECF(3)-X_j(3))/rhol_j];
         Phi1 = (obsr(p,(t-1)*5+2)-obs(p,(t-1)*5+2)...
            -obsr(p,(refprn-1)*5+2)+ ...
            obs(p,(refprn-1)*5+2))*lambda1;
         b1(tt,:) = Phi1-lambda1*Amb(tt,1)-rhok_i+rhok_j+rhol_i-rhol_j; 
      end; %t
      dx_p = A\b1;
      res_p = A*dx_p-b1;
      sigma_p = sqrt(res_p'*C*res_p/noprn);
      %Simple test for outliers
      if sigma_p < 3*sqrt(Var)
         xx = [xx dx_p]; %accumulation of solutions
         varsum =varsum+sigma_p^2;
         ressum = [ressum res_p]; %accumulation of residuals
         Normal = Normal+A'*C*A;
         RightSide = RightSide+A'*C*b1;
         t1 = t1+1;
      end % sigma_p
   end; % p
   dx = inv(Normal)*RightSide;
   Var = varsum/t1;
   S = Var*inv(Normal);
   % X_j = X_j+dx;
   sigma_X = sqrt(S(1,1));
   sigma_Y = sqrt(S(2,2));
   sigma_Z = sqrt(S(3,3));
   rho = corrcoef(S);
   
   figure((slides-1)*2+1);
   subplot(4,1,1), plot(xx(1,:)')
   switch slides 
   case 1
      title(['Estimated Baseline, \it{ L\rm_1} \rmcode only,',...
            ' float ambiguities'])
   case 2   
      title(['Estimated Baseline, \it{ L\rm_1} \rmcode only,',...
            ' fixed ambiguities'])
   case 3   
      title(['Estimated Baseline, \it{ L\rm_1} \rmcode only,',...
            ' fixed ambiguities, one changed by 1'])
   end
   ylabel('\itx \rm[m]')
   subplot(4,1,2), plot(xx(2,:)')
   ylabel('\ity \rm[m]')
   subplot(4,1,3), plot(xx(3,:)')
   ylabel('\itz \rm[m]')
   
   for i =1:size(xx,2)
      normb(i,1) = norm(xx(:,i));
   end
   subplot(4,1,4), plot(normb)
   ylabel('norm of baseline')
    
   figure((slides-1)*2+2);
   subplot(4,1,1), plot(1000*ressum(:,2:501)')
     switch slides 
   case 1
      title(['Estimated Baseline, \it{ L\rm_1} \rmcode only,',...
            ' float ambiguities'])
   case 2   
      title(['Estimated Baseline, \it{ L\rm_1} \rmcode only,',...
            ' fixed ambiguities'])
   case 3   
      title(['Estimated Baseline, \it{ L\rm_1} \rmcode only,',...
            ' fixed ambiguities, one changed by 1'])
   end
   ylabel('innovation [mm]')
   subplot(4,1,2), plot(normb)
   ylabel('norm of baseline [m]')
   subplot(4,1,3)
   plot(xcorr(1000*(ressum(4,2:501)))); 
   set(get(gca,'ylabel'),'string','corr. fnc. [mm\rm{^2}]')
   subplot(4,1,4)
   semilogy(2:998,fft(abs(xcorr(1000*ressum(4,2:500))')))
   set(get(gca,'ylabel'),'string','power spectrum [dB]')
   xlabel(['\sigma_{\rmbaseline} = ' num2str(std(1000*normb),3) ' mm'])
end
toc
%%%%%%%% end makebase.m %%%%%%%%%%%%%%%%%%%%
