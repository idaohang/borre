function makebase();
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

%Kai Borre, February 20, 1999

%THIS SAMPLE CODE DOES NOT ACCOUNT FOR CYCLE SLIPS
global Eph
% Initial computations of constants
v_light = 299792458; 	     % vacuum speed of light m/s
f1 = 154*10.23E6;		        % L1 frequency Hz
f2 = 120*10.23E6;			     % L2 frequency Hz
lambda1 = v_light/f1;	     % wavelength on L1:  .19029367  m
lambda2 = v_light/f2;	     % wavelength on L2:  .244210213 m
alpha1 = f1^2/(f1^2-f2^2);   % for use in ionospheric-free comb.
alpha2 = -f2^2/(f1^2-f2^2);

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
kprn = find(ME > 15); 
kprn = setdiff(kprn,refprn);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start of ambiguity estimation, see pages 489--490
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m,n] = size(obs(:,1));
% Definition of filter matrix, see equation (15.9)
A = [1 0       0;
     1 lambda1 0;
     1 0       0;
     1 0       lambda2];
sd = [0.3 0.005 0.3 0.005]; % standard deviation of observations
W = diag(1./sd.^2);         % diagonal weight matrix
AW = A'*W;
N = AW*A;
ef = 6;                     % first epoch
el = m;                     % last epoch
rk = 0;

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
   for i = ef:el		% i runs over epochs  
      RS = AW*b(i,:)';
      N22 = N22+N(2:3,2:3)-N(2:3,1)*N(2:3,1)'/N(1,1);
      RS21 = RS21+RS(2:3,1)-N(2:3,1)*RS(1,1)/N(1,1);
   end
   x = inv(N22)*RS21;
   K1 = round(x(1)-x(2));
   K2 = round(60*x(1)-77*x(2));
   trueN2 = round((60*K1-K2)/17);
   trueN1 = round(trueN2+K1);
   fprintf(['\nAmbiguities for PRNs'...
                 '%3.0f and %3.0f'],prns(refprn),prns(k));     
   fprintf('\nResult: N_L1: %10.0f, N_Lw: %10.0f',...
                                   trueN1, trueN1-trueN2); 
   rk = rk+1;
   amb(rk,:) = [trueN1 trueN2];
end
fprintf('\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATION OF VECTOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Computing WGS 84 coordinates of master site
Eph = get_eph('edata.dat');
time = BM(1,1)+4*epoch_interval; % we use data from 5th epoch,380600
pr = obs(5,([kprn]-1)*5+1);
pos = b_point4(pr',prns(kprn),time);
X_i = pos(1:3);
x = zeros(3,1); %Best guess for preliminary vector components
X_j = X_i;
[phi_i,lambda_i,h_i] = ...
   togeod(6378137,298.257223563,X_i(1),X_i(2),X_i(3));

% Calculation of weight matrix
noprn = length(kprn);
D = [ones(noprn,1) -eye(noprn) -ones(noprn,1) eye(noprn)];
C = inv(D*D');

% We relate each satellite with a fixed ephemeris
time = BR(1,1);
for t = 1:noprn
   col_Eph(t) = find_eph(Eph,prns(kprn(t)),time);
end
% and the reference satellite
col_Eph_r = find_eph(Eph,prns(refprn),time);

ressum = zeros(noprn,1)*inf;
Var = 100;  % Initial value for variance of an epoch solution

% Iteration for vector estimation
for iter = 1:2 % two iterations sufficient for baselines < 5 km
   [phi_j,lambda_j,h_j] =  ...
      togeod(6378137,298.257223563,X_j(1),X_j(2),X_j(3));
   Normal = zeros(3,3);
   RightSide = zeros(3,1);
   A_a = zeros(length(kprn),3);
   C_a = zeros(length(kprn),length(kprn));
   b_a = zeros(length(kprn),1);
   xx = [];
   varsum = 0;
   t1 = 0;
   firstepoch = 5; % We skip some starting epochs
   ressum = [ressum zeros(noprn,1)*inf];
   figure;
   set(gcf,'UserData',zeros(noprn,1)*inf);
   pp = plot(1,zeros(noprn,1)*inf,'.',...
      'EraseMode','None','MarkerSize',5);   
   ylabel('Residual  [m]')
   xlabel('Epochs, epoch interval 20 s')
   if iter == 1
      title('1st Iteration')
   else
      title('2nd Iteration')
   end;
   scale = 0.3;
   if iter > 1, scale = 0.03; end
   axis([0 epochs -scale scale]);
   
   for p = firstepoch:epochs
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
         % Tropospheric correction of phases, standard met. parameters
         if iter == 1
            t_cor = 0; 
         else
            [az,el_ki,d] = topocent(X_i,Xk_ECF-X_i);
            [az,el_li,d] = topocent(X_i,Xl_ECF-X_i);
            [az,el_kj,d] = topocent(X_j,Xk_ECF-X_j);
            [az,el_lj,d] = topocent(X_j,Xl_ECF-X_j);
            %el_ki, el_li, el_kj, el_lj
            t_cor = tropo(sin(el_lj*pi/180),...
               h_j*1.e-3,1013,293,50,0,0,0)...
               -tropo(sin(el_li*pi/180),....
               h_i*1.e-3,1013,293,50,0,0,0)...
               -tropo(sin(el_kj*pi/180),...
               h_j*1.e-3,1013,293,50,0,0,0)...
               +tropo(sin(el_ki*pi/180),...
               h_i*1.e-3,1013,293,50,0,0,0);
         end;
         Phi1 = (obsr(p,(t-1)*5+2)-obs(p,(t-1)*5+2)...
            -obsr(p,(refprn-1)*5+2)+ ...
            obs(p,(refprn-1)*5+2))*lambda1-t_cor;
         Phi2 = (obsr(p,(t-1)*5+4)-obs(p,(t-1)*5+4)...
            -obsr(p,(refprn-1)*5+4)+...
            obs(p,(refprn-1)*5+4))*lambda2-t_cor;
         %As observation we use the ionosphere-free combination      
         b1(tt,:) = [alpha1*(Phi1-lambda1*amb(tt,1))+ ...
               alpha2*(Phi2-lambda2*amb(tt,2))- ...
               rhok_i+rhok_j+rhol_i-rhol_j]; 
      end;
      
      dx_p = A\b1;
      res_p = A*dx_p-b1;
      sigma_p = sqrt(res_p'*C*res_p/noprn);
      % Simple test for outliers
      if (iter == 1) | (sigma_p < 3*sqrt(Var)) 
         xx = [xx dx_p];
         varsum = varsum+sigma_p^2;
         ressum = [ressum res_p];
         Normal = Normal + A'*C*A; 
         RightSide = RightSide + A'*C*b1; 
         A_a = A_a+A;
         C_a = C_a+C;
         b_a = b_a+b1;
         t1 = t1+1;
      end;
      set(pp,'XData',p*ones(noprn,1),'YData',res_p) 
      drawnow
   end;  % end p
   dx = inv(Normal)*RightSide;
   Var = varsum/t1;
   S = Var*inv(Normal);
   X_j = X_j+dx;
   sigma_X = sqrt(S(1,1));
   sigma_Y = sqrt(S(2,2));
   sigma_Z = sqrt(S(3,3));
   fprintf('\n');
   fprintf('Iteration # %2.0f\n', iter);
   fprintf(['Correction x, y, z: %8.4f'....
                         ' %8.4f %8.4f\n'], dx(1),dx(2),dx(3));
   fprintf(['Sigma x, y, z:      %8.4f'....
                   ' %8.4f %8.4f\n'], sigma_X,sigma_Y,sigma_Z);
   fprintf('Epochs rejected %2.0f\n', epochs-firstepoch+1-t1);
   rho = corrcoef(S);
end; % end iter

% Correction for antenna heights
fids = fopen('sdata.dat');
hmr = fread(fids,2,'double');
%master - rover
[up(1,1),up(2,1),up(3,1)] = enu2xyz(phi_i,lambda_i,0,0,hmr(2)-hmr(1));
vect = X_j-X_i + up;
fprintf(['\nFINAL VALUES FOR VECTOR\n deltaX: %10.3f'....
                     ' deltaY: %10.3f deltaZ: %10.3f\n'],...
                                      vect(1),vect(2),vect(3));
res_a =A_a*dx-b_a;
Sigma = norm(res_a)^2*inv(A_a'*C_a*A_a);   
fprintf('Sigma_x  %10.3f Sigma_y %10.3f Sigma_z %10.3f\n\n',...
            sqrt(Sigma(1,1)),sqrt(Sigma(2,2)),sqrt(Sigma(3,3)));
fprintf('\n Distance %10.3f\n', norm(vect));

figure(3);
plot(xx'*1.e3)
ylabel('Residuals in x, y, z [mm]')
title('Vector Estimation');

figure(4);
plot(ressum(1,:)*1.e3)
ylabel(['L1 Residuals, Sv ',...
      num2str(prns(refprn)),' - Sv ',num2str(prns(1)),' [mm]'])
title(['Vector Estimation Through All ',num2str(iter),' Iterations']);

figure(5);
plot(ressum(noprn,:)*1.e3)
ylabel(['L2 Residuals, Sv ',...
      num2str(prns(refprn)),' - Sv ',num2str(prns(1)),' [mm]'])
title(['Vector Estimation Through All ',num2str(iter),' Iterations']);

figure(6);  
plot(ressum')
title('Baseline Estimation From DD Code & Phase Obsv.s')
ylabel('Residuals [m]')
xlabel(['Epochs, ',num2str(iter),' Iterations'])

fprintf('\nElapsed time (sec): %3.2f\n', toc);

%----------------------------------------------

function [x, y, z] = enu2xyz(phi,lambda,e,n,u);
%ENU2XYZ  Transformation of [e;n;u] vector from local to geocentric
%   	    system. The local system has origin at (phi,lambda)

dtr = pi/180;
cl = cos(lambda*dtr);  sl = sin(lambda*dtr);
cb = cos(phi*dtr);	  sb = sin(phi*dtr);
F = [-sl -sb*cl cb*cl;
      cl -sb*sl cb*sl;
      0	  cb      sb];
global_vector = F*[e; n; u];
x = global_vector(1);
y = global_vector(2);
z = global_vector(3);

%%%%%%%%% end makebase.m %%%%%%%%%%%%%%%%%%%%
