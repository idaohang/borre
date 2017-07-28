%function makemult()
%MAKEMULT  We assume that you already have run the M-file bdata 
%          with b-files from master and rover. The call would be
%             bdata('masterfile','roverfile') or
%	           bdata('b0810a94.076','b0005a94.076')
%          That creates the file bdata.dat.
%          Likewise we assume that the calls
%             edata('masterfile','roverfile') and
%             sdata('masterfile','roverfile')
%          have created edata.dat and sdata.dat.
%        
%          Estimating multipath for L1 and L2 code.
%
%          References refer to
%          Strang, Gilbert and Borre, Kai (19997): Linear Algebra,
%                     Geodesy, and GPS. Wellesley-Cambridge Press.

%Kai Borre, February 17, 1999

%THIS SAMPLE CODE DOES NOT ACCOUNT FOR CYCLE SLIPS

% Initial computations of constants
v_light = 299792458; 	     % vacuum speed of light m/s
f1 = 154*10.23E6;		        % L1 frequency Hz
f2 = 120*10.23E6;			     % L2 frequency Hz
lambda1 = v_light/f1;	     % wavelength on L1:  .19029367  m
lambda2 = v_light/f2;	     % wavelength on L2:  .244210213 m
beta = (f1/f2)^2;

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
kprn = find(ME > 15); 
kprn = setdiff(kprn,refprn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multipath estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The observation vector contains double differences 
%of     P1   Phi1   P2    Phi2. All data in a row pertain to
%the same epoch
%Master
P1_i = repair_c(obs(:,(refprn-1)*5+1));
M1_ik = P1_i-(beta+1)/(beta-1)*obs(:,(refprn-1)*5+2)*lambda1+...
                        2/(beta-1)*obs(:,(refprn-1)*5+4)*lambda2;
M1_ik = M1_ik-M1_ik(1,1);
M1_run = [];
for k = kprn  
    P1_l = repair_c(obs(:,(k-1)*5+1));
    M1_il = P1_l-(beta+1)/(beta-1)*obs(:,(k-1)*5+2)*lambda1+ ...
                            2/(beta-1)*obs(:,(k-1)*5+4)*lambda2;
    M1_il = M1_il-M1_il(1,1);
    M1_run = [M1_run M1_il];
end


%Rover
beg = 3; % first epoch usable
P1_j = repair_c(obsr(beg:epochs,(refprn-1)*5+1));
M1_jk = P1_j-(beta+1)/(beta-1)*obsr(beg:epochs,(refprn-1)*5+2)*...
         lambda1+2/(beta-1)*obsr(beg:epochs,(refprn-1)*5+4)*lambda2;
M1_jk = M1_jk-M1_jk(1,1);
M1_runr = [];

for k = kprn  
    P1_jl = repair_c(obsr(beg:epochs,(k-1)*5+1));
    M1_jl = P1_jl-(beta+1)/(beta-1)*obsr(beg:epochs,(k-1)*5+2)*...
               lambda1+2/(beta-1)*obsr(beg:epochs,(k-1)*5+4)*lambda2;
    M1_jl = M1_jl-M1_jl(1,1);
    M1_runr = [M1_runr M1_jl];
end

plot(M1_ik,'Linewidth',2)
hold on
plot(M1_run)
title('P-code multipath on L1 from base receiver')
ylabel('P-code multipath [m]')
xlabel('Epochs, epoch interval 20 s')

%Single differences of multipath on L1 at master
figure;
plot(M1_run-M1_ik*ones(1,5))
title(['P-code multipath on L1, all PRN''s minus ref. PRN, '...
              'base receiver'])
ylabel('P-code multipath [m]')
xlabel('Epochs, epoch interval 20 s')

%Single differences of multipath on L1 at rover
figure;
plot(M1_runr-M1_jk*ones(1,5))
title(['P-code multipath on L1, all PRN''s minus ref. PRN,'...
               'rover receiver'])
ylabel('P-code multipath [m]')
xlabel('Epochs, epoch interval 20 s')


%Double differences of multipath on L1
figure;
plot(M1_run(beg:epochs,:)-M1_runr-...
               (M1_ik(beg:epochs)-M1_jk)*ones(1,5))
title('P-code multipath on L1, Double differences')
ylabel('P-code DD multipath [m]')
xlabel('Epochs, epoch interval 20 s')

%%%%%%%%%%%%%%%%%%%% end makemult.m  %%%%%%%%%%%%%%%%%%