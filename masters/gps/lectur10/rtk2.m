%RTK    Script for developing RTK code

% Copyright by Kai Borre
% March 18, 2003

global Eph ambi ty
global ATA ATb
global X_ista1 X_ista2
global sat_index

% Initial computations of constants
v_light = 299792458; 	     % vacuum speed of light m/s
f1 = 154*10.23E6;		     % L1 frequency Hz
f2 = 120*10.23E6;		     % L2 frequency Hz
lambda1 = v_light/f1;	     % wavelength on L1:  .19029367  m
lambda2 = v_light/f2;	     % wavelength on L2:  .244210213 m
beta = (f1/f2)^2;

ty = 30;       % Maximum number of normal equations 
ATA = zeros(ty,ty);
ATA = ATA+1.e-20*eye(ty);
ATA(4,4) = 1.e16; 
ATb = zeros(ty,1);
sat_index = zeros(25,1); % is used in locate.m
ambi = zeros(30,1); 

% Getting the ephemerides
% In this code we treat the ephemerides the usual way.
% In case of true RTK we start by dowmloading possible
% ephemerides via the COM port. Any new ephemeris knocks
% on the door itself. Thus it becomes easy, without too
% much bookkeeping to store only the lastest ephemerides
% for all tracked saltellies.
Eph = edata('ELT03A97.149');

% Index 1 relates to Master receiver, 2 to Rover receiver
[rawheader1,fidobs1] = z12head('BLT03A97.149');
rawnav1 = z12nav(fidobs1);
time1 = rawnav1.rcv_time;

[rawheader2,fidobs2] = z12head('BLT04A97.149');
rawnav2 = z12nav(fidobs2);
time2 = rawnav2.rcv_time;

% The following matching of starting times for both nav-files is only
% needed when we are reading from files--when data arrive in real time
% they are expected to be simultaneous
if (time2 < time1),
    while 1
        rawnav2 = z12nav(fidobs2);
        time2 = rawnav2.rcv_time;
        if time2 == time1, break, end;
    end
else
    while 1
        rawnav1 = z12nav(fidobs1);
        time1 = rawnav1.rcv_time;
        if time1 == time2, break, end;
    end
end;

%Correction for different antenna heights
% Rover = ant_delta(1), Master = ant_delta(2)
ant_delta = sdata('slt03a97.149','slt04a97.149');

rawnav1 = z12nav(fidobs1);
time1 = rawnav1.rcv_time;
rawnav2 = z12nav(fidobs2);
time2 = rawnav2.rcv_time;
sv1 = rawnav1.nest(1).svprn;
% Deleting rows for non-tracked PRNs
sv1 = sv1(~isnan(sv1));
el1 =rawnav1.nest(1).elevation([sv1]);

% Delete Sv with elevation smaller than 10 degrees
sv1(el1 < 10) = []; 
sv2 = rawnav2.nest(1).svprn;
sv2 = sv2(~isnan(sv2));
sv = intersect(sv1,sv2);

% We read and skip some more epochs 
for p = 1:20
    rawnav1 = z12nav(fidobs1);
    rawnav2 = z12nav(fidobs2);
end       

% Get P1 code (nest(2)) observations from Master (1) and Rover (2)
obs1 = rawnav1.nest(2).rawrange([sv],1);
obs2 = rawnav2.nest(2).rawrange([sv],1);
% We compute the Master and Rover positions from P1 code
[X1, el] = recpo_ls(obs1, sv, time1, Eph);
X_ista1 = X1(1:3,1);
[X2, el] = recpo_ls(obs2, sv, time2, Eph);
X_ista2 = X2(1:3,1);
%[phi_i,lambda_i,h_i] =  ...
%    togeod(6378137,298.257223563,X_ista1(1),X_ista1(2),X_ista1(3));

base = X_ista2-X_ista1

for p = 1:50
    sv1 = rawnav1.nest(1).svprn;
    sv1 = sv1(~isnan(sv1));
    el1 =rawnav1.nest(1).elevation([sv1]);
    % Delete low PRNs
    sv1(el1 < 10) = []; 
    sv2 = rawnav2.nest(1).svprn;
    sv2 = sv2(~isnan(sv2));
    sv = intersect(sv1,sv2);
    noprn = length(sv);
    
    rawnav1 = z12nav(fidobs1);
    time1 = rawnav1.rcv_time;
    rawnav2 = z12nav(fidobs2);
    time2 = rawnav2.rcv_time;
    % We continue to check for missing epochs of data
    while time1 ~= time2
        if time1 > time2
            rawnav2 = z12nav(fidobs2);
            time2 = rawnav2.rcv_time;
            disp('read 2')
        end
        if time2 > time1
            rawnav1 = z12nav(fidobs1);
            time1 = rawnav1.rcv_time;
            disp('read 1')
        end
        if time1 == time2, break, end
    end
    % Rawrange on P1 at Master 1 and Rover 2 is now in
    % obs1 and obs2, sorted corresponding to the 
    % sequence of PRNs in sv
    obs1p1 = rawnav1.nest(2).rawrange([sv],1); 
    obs1p2 = rawnav1.nest(3).rawrange([sv],1);  
    obs2p1 = rawnav2.nest(2).rawrange([sv],1);     
    obs2p2 = rawnav2.nest(3).rawrange([sv],1);     
    % Next we read the phases for Master on L1 and L2 from nest(2) and nest(3)
    obs1ph1 = rawnav1.nest(2).carphase([sv],1);  % Phi1 at Master
    obs1ph2 = rawnav1.nest(3).carphase([sv],1);  % Phi2 at Master
    obs2ph1 = rawnav2.nest(2).carphase([sv],1);  % Phi1 at Rover
    obs2ph2 = rawnav2.nest(3).carphase([sv],1);  % Phi2 at Rover
    sum_norm(time1,sv,[obs1ph1,obs1ph2,obs1p1,obs1p2],...
        [ant_delta(2), 0, 0], time2,sv,[obs2ph1,obs2ph2,obs2p1,obs2p2], ...
        [ant_delta(1),0,0],lambda1)        
    
end
    % We delete the zero part of the normals and find the solution
    cut = max(find(diag(ATA > 0.0000001)));
    ATA(cut+1:20,:) = [];
    ATA(:,cut+1:20) = [];
    ATb(cut+1:20) = [];
    M = inv(ATA);
    % The recent Lambda implementation is very sensitive to non-symmetry,
    % so we repair
    M = (M+M')/2;
    dSol = M*ATb;
    trace = M(1,1)+M(2,2)+M(3,3);
    pdop = sqrt(trace);
    
    %LAMBDA decorrelation: We compute the integer ambiguities N,
    %       and correct the previous solution
    fidlog = fopen('lambda.log','wt');
    [N,disall,Qcheck,Z] = lambda1(dSol(5:cut,1),M(5:cut,5:cut),cut-4,fidlog);
    % We choose the first column of N, out the (six) possible.
    % We may compute the success rate by calling: success.m 
    dSol = inv(ATA(1:3,1:3))*(ATb(1:3,1)-ATA(1:3,5:cut)*N(:,1));
    baseline = X_ista2-X_ista1-dSol(1:3,1);
    fprintf('\n Vector from 1 to 2:%12.3f  %12.3f  %12.3f\n',...
        baseline(1), baseline(2), baseline(3)) 
    
fclose('all');
%%%%%%%%%%%%%%%%%%%%%%%% end rtk.m  %%%%%%%%%%%%%%
