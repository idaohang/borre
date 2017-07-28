function[Pos,amb,goodm,goodr,P1goodm,iprnm,iprnr,ref_ind] = teunissen(X_i,X_j,no_s,no_i)
%TEUNISSEN        Ambiguity estimation according to Peter Teunissen's method,
%

% Written by Kai Borre
% November 9, 2007

global EPH

% Initial computations of constants
v_light = 299792458; 	     % vacuum speed of light m/s
f1 = 154*10.23E6;		     % L1 frequency Hz
f2 = 120*10.23E6;		     % L2 frequency Hz
lambda1 = v_light/f1;	     % wavelength on L1:  .19029367  m
lambda2 = v_light/f2;	     % wavelength on L2:  .244210213 m

Pos = [];
x = zeros(3+2*(no_s-1),1);

for i = 1:no_i
    [towm,prnm,P1m,Phi1m,P2m,Phi2m] = readGrilM; % fidm
    [towr,prnr,P1r,Phi1r,P2r,Phi2r] = readGrilR; % fidr
    % Find the total number of ephemerides
    prnE = find(EPH(1,:) > 0);
    % Find number lgm of master PRNs with known ephemerides
    [goodm, iprnm] = intersect(prnm,prnE);
    lgm = length(goodm);
    % Find number lgr of rover PRNs with known ephemerides
    [goodr, iprnr] = intersect(prnr,prnE);
    % Note that the above operations leave the PRNs sorted
    % in increasing order

    % Compute elevation angle el for deleting low satellites
    [pos, el] = recposRTK(towm,prnm,P1m');
    Pos = [Pos pos];

    % Delete PRNs with elevation below 10 degrees
    goodm(el < 10) = [];
    goodr = goodm;
    % Number of deleted common PRNs due to low elevation
    No_delPRN = lgm-length(goodm);
    if i == 1
        disp(['Number of PRNs deleted due to low elevation ' num2str(No_delPRN)])
    end
    lgm = length(goodm);
    % The PRN with largest elevation is selected as reference
    [y,ref_ind] = max(el);

    % ref_ind keeps the index of the reference PRN and the
    % rest are in sat_ind
    sat_ind = setdiff(1:lgm,ref_ind);
    if lgm >= 4,
        P1goodm = [];
        Phi1goodm = [];
        P2goodm = [];
        Phi2goodm = [];
        for m = 1:lgm
            % rearranging all useful master observations
            P1goodm = [P1goodm; P1m(iprnm(m))];
            Phi1goodm = [Phi1goodm; Phi1m(iprnm(m))];
            P2goodm = [P2goodm; P2m(iprnm(m))];
            Phi2goodm = [Phi2goodm; Phi2m(iprnm(m))];
        end
        P1goodr = [];
        Phi1goodr = [];
        P2goodr = [];
        Phi2goodr = [];
        for m = 1:lgm
            % rearranging all useful rover observations
            P1goodr = [P1goodr; P1r(iprnr(m))];
            Phi1goodr = [Phi1goodr; Phi1r(iprnr(m))];
            P2goodr = [P2goodr; P2r(iprnr(m))];
            Phi2goodr = [Phi2goodr; Phi2r(iprnr(m))];
        end
    end % lgm

    % We loop over all good PRNs (sat_ind).
    % The outer i-loop runs over the epochs from which we
    % want to estimate the ambiguities; the right sides are
    % added and averaged. This works because
    % of the additive character of normal equations

    ref_PRN = goodm(ref_ind);
    tt = 0; %counts no. of PRNs
    A1 = []; 
    lgm = length(sat_ind);
    
    % Select data for the ref satellite at master
    [rhok_i,Xk_ECF] = get_rho(towm, P1goodm(ref_ind), ref_PRN, X_i(1:3));
    % Select data for the ref satellite at rover
    [rhok_j,Xk_ECF] = get_rho(towr, P1goodr(ref_ind), ref_PRN, X_j(1:3));
    for t = sat_ind
        tt = tt+1;
        [rhol_i,Xl_ECF] = get_rho(towm, P1goodm(t), goodm(t), X_i(1:3));
        [rhol_j,Xl_ECF] = get_rho(towr, P1goodr(t), goodr(t), X_j(1:3));
        A0 = ...
            [(Xk_ECF(1)-X_j(1))/rhok_j - (Xl_ECF(1)-X_j(1))/rhol_j, ...
            (Xk_ECF(2)-X_j(2))/rhok_j - (Xl_ECF(2)-X_j(2))/rhol_j, ...
            (Xk_ECF(3)-X_j(3))/rhok_j - (Xl_ECF(3)-X_j(3))/rhol_j];
        A1 = [A1; A0];
        Phi1 = Phi1goodr(t)-Phi1goodm(t)-Phi1goodr(ref_ind)+Phi1goodm(ref_ind);
        Phi2 = Phi2goodr(t)-Phi2goodm(t)-Phi2goodr(ref_ind)+Phi2goodm(ref_ind);
        b(tt,:) = Phi1-lambda1*x(3+tt,1);
        b(lgm+tt,:) = Phi2-lambda2*x(3+lgm+tt,1);
        bk(tt,:) = rhok_i-rhok_j-rhol_i+rhol_j;
        bk(lgm+tt,:) = rhok_i-rhok_j-rhol_i+rhol_j;
    end; %t, i.e. loop over PRNs of a single epoch

    N = zeros(3+2*lgm,3+2*lgm);
    rs = zeros(3+2*lgm,1);
    D = [ones(lgm,1) -eye(lgm) -ones(lgm,1) eye(lgm)];
    Sigma = D*D';
    A_modi = eye(lgm);
    A_modi(:,ref_ind) = -ones(lgm,1);
    A_aug = [A1 lambda1*A_modi 0*eye(lgm); A1 0*eye(lgm) lambda2*A_modi];
    N = N + A_aug'*kron(eye(2), Sigma)*A_aug;
    rs = rs + A_aug'*kron(eye(2),Sigma)*(b-bk);
end; % no_i

PP = pinv(N);
x = PP*rs;
[amb,sqnorm,Sigma_afixed, Z] = lambda(x(4:4+2*lgm-1,1), PP(4:4+2*lgm-1,4:4+2*lgm-1));

%%%%%%%%%%%%%%%% end teunissen.m  %%%%%%%%%%%%%%%%%%%%%%%