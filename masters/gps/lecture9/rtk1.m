% MAIN3  RTK computation of baseline using dobble differenced
%        carrier phase observations on L1.
%        Data read "as if" in real time from Topcon receivers
%        according to the GRIL format
%
% The code does not handle
%	     1. cycle slips, and
%	     2. outliers.

% Written by Kai Borre
% November 12, 2007

global EPH fidm fidr

% Initial computations of constants
v_light = 299792458; 	     % vacuum speed of light m/s
f1 = 154*10.23E6;		     % L1 frequency Hz
f2 = 120*10.23E6;		     % L2 frequency Hz
lambda1 = v_light/f1;	     % wavelength on L1:  .19029367  m
lambda2 = v_light/f2;	     % wavelength on L2:  .244210213 m
alpha = (f1/f2)^2;

EPH = zeros(21,32); %global variable to hold ephemerides
masterfile = '19jan07m.log';
roverfile =  '19jan07r.log';
fidm = fopen(masterfile);
fidr = fopen(roverfile);
eof = 0;

% Read data from one epoch to compare epoch times in the two files.
% The actual master and rover files have identical epochs.
% This step is needed because we read from files. In real
% RTK all epochs should automatically be synchronized
[towm,prnm,P1m,Phi1m,P2m,Phi2m] = readGrilM; % fidm
[towr,prnr,P1r,Phi1r,P2r,Phi2r] = readGrilR; % fidr
if towm < towr % only this loop tested
    dt = towr-towm;
    for i = 1:dt
        [towm,prnm,P1m,Phi1m,P2m,Phi2m] = readGrilM;
    end
elseif towm > towr
    dt = towm-towr;
    for i = 1:dt
        [towr,prnr,P1r,Phi1r,P2r,Phi2r] =  readGrilR;
    end
elseif (round(towm) == round(towr))
    [towm,prnm,P1m,Phi1m,P2m,Phi2m] = readGrilM;
    [towr,prnr,P1r,Phi1r,P2r,Phi2r] = readGrilR;
end

%Computing WGS 84 coordinates of master site
[X_i, El] = recposRTK(towm,prnm,P1m');
[phi_i,lambda_i,h_i] = ...
    togeod(6378137,298.257223563,X_i(1),X_i(2),X_i(3));
fprintf('\nStart of Common Epoch Time:  %8.1f  %8.1f',towm,towr)
fprintf('\n')

% We read no_i epoch to estimate the ambiguities for L1 and
% L2 phase observations. The estimate is based on a method
% indicated by Clyde C. Goad

no_i = 30;
%%[Pos,amb,goodm,goodr,P1goodm,iprnm,iprnr,ref_ind] = goad(no_i);
%%ref_PRN = goodm(ref_ind);
[Pos,amb,goodm,goodr,P1goodm,iprnm,iprnr,ref_ind] = teunissen(X_i, X_i, length(prnm), no_i);
ref_PRN = goodm(ref_ind);
lgm = length(goodm);
D = [ones(lgm-1,1) -eye(lgm-1) -ones(lgm-1,1) eye(lgm-1)];
Sigma = D*D';

figure(1);
plot((Pos(1:3,:)-Pos(1:3,1)*ones(1,no_i))')
title(['Variation of Master Station Coordinates for First ' ...
    num2str(no_i) ' Epochs [m]'])

%Initially, the rover site equals the master site
X_j = X_i;
phi_j = phi_i;
lambda_j = lambda_i;
h_j = h_i;

%Preparations for Figure 2
figure(2);
set(gcf,'UserData',zeros(3,1)*inf);
pp = plot(1,zeros(3,1)*inf,'.',...
    'Erasemode','None','MarkerSize',2);
set(pp(1),'Color','black','Marker','*')
set(pp(2),'Color','red','Marker','o')
set(pp(3),'Color','blue')
ylabel('Innovation in Vector Components [m]')
xlabel('Epochs, interval 1 s')
scale = 1; %%0.05;
epochs = 2550; % A priori knowledge
axis([0 epochs -scale scale]);
figdata = [];

%Preparation for position filter
P = 10^8*eye(4);        % covariance for state vector
Q = 0.05^2*eye(4);	     % covariances of system
R = 0.005^2*inv(Sigma);	 % covariances of observations
A = [];
p = 0; % number of epochs processed
x = zeros(4,1);

%Processing epoch-by-epoch
while 1
    p = p+1;
    % We read an epoch of master data until end of file
    [towm,prnm,P1m,Phi1m,P2m,Phi2m,eof] = readGrilM;
    if eof == 1, break, end
    % In the obervation file any missing observation was
    % substituted by an NaN; now we test if NaN is present and
    % delete the corresponding PRN and pertinent data
    mis_obsm = isnan(P1m);
    if sum(mis_obsm) > 0
        fprintf('\nDeleted Master Observation in Epoch %4.0f ',towm);
        colm = find(mis_obsm == 1);
        prnm(:,colm) = [];
        P1m(:,colm) = [];
        Phi1m(:,colm) = [];
        P2m(:,colm) = [];
        Phi2m(:,colm) = [];
    end

    % We read an epoch of rover data
    [towr,prnr,P1r,Phi1r,P2r,Phi2r] = readGrilR;
    % We test if NaN is present and delete the corresponding
    % prn and its data
    mis_obsr = isnan(P1r);
    if sum(mis_obsr) > 0
        fprintf('\nDeleted Rover Observation in Epoch %4.0f ',towr);
        colr = find(mis_obsr == 1);
        prnr(:,colr) = [];
        P1r(:,colr) = [];
        Phi1r(:,colr) = [];
        P2r(:,colr) = [];
        Phi2r(:,colr) = [];
    end

    % A new PRN most often rises at different epochs at master
    % and rover. Hence we need to delete the observations from
    % one receiver until the PRN is given with complete data
    % from both receivers
    if length(prnm) ~= length(prnr)
        [prnm_cut,im,ir] = intersect(prnm,prnr);
        [new_prn,new_i] = setdiff(prnm,prnm_cut);

        % Deleting additional PRN data from master data
        if length(prnm) > length(prnr)
            fprintf('\n Deleted master data from a PRN')
            prnm(:,new_i) = [];
            P1m(:,new_i) = [];
            Phi1m(:,new_i) = [];
            P2m(:,new_i) = [];
            Phi2m(:,new_i) = [];
        end
        % Deleting additional PRN data from rover data
        if length(prnr) > length(prnm_cut)
            fprintf('\n Deleted rover data from a PRN')
            prnr(:,new_i) = [];
            P1r(:,new_i) = [];
            Phi1r(:,new_i) = [];
            P2r(:,new_i) = [];
            Phi2r(:,new_i) = [];
        end
    end

    if length(goodm) < 4
        fprintf('\nNot sufficient number of PRNs in epoch %8.0 ',towm,towr)
        continue
    end

    % prnm is a list of PRNs at master receiver with known
    % ephemeris and complete data set
    lprnm = length(prnm);
    % prnr is a list of PRNs at rover receiver with known
    % ephemeris and complete data set
    lprnr = length(prnr);

    if lprnm >= 4
        if lprnm > lgm % lgm holds the number of PRNs in the previous epoch
            if lprnm ~= lprnr, continue, end
            fprintf('\nWe estimate AMB anew, at epoch %6.0f\n', towm)
            old_amb = amb;
            [Pos,amb,goodm,goodr,P1goodm,iprnm,iprnr,ref_ind] = goad(no_i);
            ref_PRN = goodm(ref_ind);   
          %%  [Pos,amb,goodm,goodr,P1goodm,iprnm,iprnr] = teunissen(X_i,X_j,lprnm,no_i);
            ref_PRN = goodm(ref_ind);
            new_amb = amb(size(old_amb,1)+1:end,:);
            amb = [old_amb; new_amb];
            lgm = length(goodm);
            D = [ones(lgm-1,1) -eye(lgm-1) -ones(lgm-1,1) eye(lgm-1)];
            Sigma = D*D';
            R = 0.005^2*inv(Sigma);	 % covariances of observations
        end
        % re-ordering all useful master observations in
        % increasing order
        P1goodm = [];
        Phi1goodm = [];
        P2goodm = [];
        Phi2goodm = [];
        for m = 1:lprnm
            P1goodm = [P1goodm; P1m(iprnm(m))];
            Phi1goodm = [Phi1goodm; Phi1m(iprnm(m))];
            P2goodm = [P2goodm; P2m(iprnm(m))];
            Phi2goodm = [Phi2goodm; Phi2m(iprnm(m))];
        end
        % re-ordering all useful rover observations
        % to match the master sequence
        P1goodr = [];
        Phi1goodr = [];
        P2goodr = [];
        Phi2goodr = [];
        for m = 1:lprnr
            P1goodr = [P1goodr; P1r(iprnr(m))];
            Phi1goodr = [Phi1goodr; Phi1r(iprnr(m))];
            P2goodr = [P2goodr; P2r(iprnr(m))];
            Phi2goodr = [Phi2goodr; Phi2r(iprnr(m))];
        end
    end % lgm

    sat_ind = 1:lgm;
    sat_ind(ref_ind) = [];
    % Select data for the ref satellite at master
    [rhok_i,Xk_ECF] = get_rho(towm, P1goodm(ref_ind), ref_PRN, X_i(1:3));
    % Select data for the ref satellite at rover
    [rhok_j,Xk_ECF] = get_rho(towr, P1goodr(ref_ind), ref_PRN, X_j(1:3));

    tt = 0; %counts no. of PRNs
    for t = sat_ind
        tt = tt+1;
        [rhol_i,Xl_ECF] = get_rho(towm, P1goodm(t), goodm(t), X_i(1:3));
        [rhol_j,Xl_ECF] = get_rho(towr, P1goodr(t), goodr(t), X_j(1:3));
        A(tt,:) = ...
            [(Xk_ECF(1)-X_j(1))/rhok_j - (Xl_ECF(1)-X_j(1))/rhol_j, ...
            (Xk_ECF(2)-X_j(2))/rhok_j - (Xl_ECF(2)-X_j(2))/rhol_j, ...
            (Xk_ECF(3)-X_j(3))/rhok_j - (Xl_ECF(3)-X_j(3))/rhol_j, 1];
        % Tropospheric correction of phases, standard met. parameters
        [az,el_ki] = topocent(X_i(1:3),Xk_ECF-X_i(1:3));
        [az,el_li] = topocent(X_i(1:3),Xl_ECF-X_i(1:3));
        [az,el_kj] = topocent(X_j(1:3),Xk_ECF-X_j(1:3));
        [az,el_lj] = topocent(X_j(1:3),Xl_ECF-X_j(1:3));
        t_cor = tropo(sin(el_lj*pi/180),...
            h_j*1.e-3,1013,293,50,0,0,0)...
            -tropo(sin(el_li*pi/180),....
            h_i*1.e-3,1013,293,50,0,0,0)...
            -tropo(sin(el_kj*pi/180),...
            h_j*1.e-3,1013,293,50,0,0,0)...
            +tropo(sin(el_ki*pi/180),...
            h_i*1.e-3,1013,293,50,0,0,0);
        Phi1 = Phi1goodr(t)-Phi1goodm(t)-Phi1goodr(ref_ind)+Phi1goodm(ref_ind)-t_cor;
        Phi2 = Phi2goodr(t)-Phi2goodm(t)-Phi2goodr(ref_ind)+Phi2goodm(ref_ind)-t_cor;
        I = Phi2-lambda2*amb(tt,2)-Phi1+lambda1*amb(tt,1);
        b(tt,:) = Phi1-rhok_i+rhok_j+rhol_i-rhol_j+lambda1*amb(tt,1)+I; 
    end; %t, i.e. loop over PRNs of a single epoch

    %Extended filter, cf. page 509--510
    P = P+Q;
    K = P*A'*inv(A*P*A'+R);
    x = K*b;
    X_j = X_j+x;
    P = (eye(size(x,1))-K*A)*P;
    set(pp,'XData',p*ones(3,1),'YData',x(1:3))
    drawnow
    figdata = [figdata X_j(1:3,1)];
    [phi_j,lambda_j,h_j] = ...
        togeod(6378137,298.257223563,X_j(1),X_j(2),X_j(3));
end; %  p, i.e. loop over epochs

xf = X_j-X_i; %+up;
fprintf(['\n\nFINAL VALUES FOR VECTOR\n dX: %8.3f'....
    ' dY: %8.3f dZ: %8.3f\n'],...
    xf(1),xf(2),xf(3));
fprintf('\nDistance %10.3f\n', norm(xf(1:3)));

figure(3);
plot((figdata(:,1)*ones(1,size(figdata,2))-figdata)')
ylabel('[m]')
legend('\itX','\itY','\itZ')
%%%%%%%%%%%%%%%%%%%%%% end main3.m  %%%%%%%%%%%%%%%%%%%