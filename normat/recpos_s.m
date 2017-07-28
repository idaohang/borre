function [phi,lambda] = recpos_s(B,tcorr)
% RECPOS_S      Least-squares searching for receiver position. Given 
%               4 or more pseudoranges and data from ephemerides. 
%               Idea to this script originates from Clyde C. Goad
%
% Made by Kai Borre 04-19-96
% Copyright (c) by Kai Borre
% Revision 1.0 Date: 1997/09/24
% Modified by group 923 AAU, Date: 1998/11/25
% Improved by: 1) A better and faster intial guess taking 
%                 care of the problem near the Greenwich meridian
%              2) A simpler data input format
%              3) Using rotation outside loops
%              4) Using a faster/converging grid + an optimized
%                 interation bound
%              5) Removing all sin/cos ops. from code, only those
%                 which can be tabelized remain
%              6) Moving asin ops. out of loops leaving only two
%                 asins left 

XS = B(:,1:3)';
pseudoranges = B(:,4);
tcorr = zeros(length(pseudoranges)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define constants and initialize variables % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grid parameters
iterations = 38;     % iterations in main loop, signifies precision
scalediv = 1.5;      % scale division factor  
ndiv = 1;            % no. of slices in inclination of spherical cap 
azimuth_angle = 60;  % angle to slice spherical cap in azimuth 
last_angle = 300;    % maximum angle in azimuth
scale = 70;          % initialize scale (deg from zenith to cap edge) 

% Universal and global constants 
dtr = pi/180;                     % conversion factor between degrees and radians
vlight = 299792458;	              % vacuum speed of light in m/s
Omegae_dot = 7.292115147e-5;      % Earth rotation rate in rad/s
semi_major_axis = 6378137;        % Earth semi-major axis       
flattening = 1/298.257223563;     % Earth flattening
esq = (2-flattening)*flattening;  % Earth eccentricity squared (6.69437999014e-3)

% Inverse constants, used to avoid division
inv_scalediv = 1/scalediv;     % inverse of scalediv (used only with psi LUT's) 
inv_ndiv = 1/ndiv;             % inverse of ndiv     (used only with psi LUT's) 
inv_vlight = 1/vlight;         % inverse of vlight     
inv_dtr = 180/pi;              % inverse of dtr


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate lookup tables for sine and cosine of alpha and psi, respectively %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_counter = 0;
for alpha = 0:azimuth_angle:last_angle     % generate cos/sin table for alpha (2*6 values)
    alpha_counter = alpha_counter+1;       % inc alpha_counter
    cos_alpha(alpha_counter) = cos(alpha*dtr);      
    sin_alpha(alpha_counter) = sin(alpha*dtr);      
end

psi_counter = 0;
for iter = 1:iterations
    for b = 1:ndiv                           % generate cos/sin table for psi (2*38 values)
        psi_counter = psi_counter+1;         % inc psi_counter
        psi = b*scale*inv_ndiv;              % angle to zenith in spherical cap
        cos_psi(psi_counter) = cos(psi*dtr); % compute sine and cosine of psi
        sin_psi(psi_counter) = sin(psi*dtr); %
        scale = scale*inv_scalediv;          % reduce scale
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute a start coordinate to get the algorithm going % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Space vehicle positions are known in the ECEF system in form of 
% Cartesian coordinates (X,Y,Z). The pseudoranges and tcorr are 
% known as well.

% Initialize phi_old & lambda_old thus establishing the initial 
% receiver position estimate using the known Cartesian coordinates
% (phi_old,lambda_old) signifies the initial point (zenith) of 
% the present spherical cap throughout the algorithm

% compute mass mean of the SV positions and find the appropriate ECEF octant
signs = sign(mean(XS'));            % define a vector of signums
signs(find(signs)==0) = 1;          % repair sign if input==0
% cos(phi) is constantly sqrt(1/2) and sin(phi) has same sign as Z
% cos(lambda) has same sign as X and sin(lambda) has same sign as Y 
cos_phi_old = sqrt(1/2);
if signs(3)==1, sin_phi_old = sqrt(1/2);    else sin_phi_old = -sqrt(1/2);    end;
if signs(1)==1, cos_lambda_old = sqrt(1/2); else cos_lambda_old = -sqrt(1/2); end;
if signs(2)==1, sin_lambda_old = sqrt(1/2); else sin_lambda_old = -sqrt(1/2); end;

% Convert sine and cosine of initial position estimate to 
% Cartesian ECEFs. Beginning of simplified FRGEOD.M unfolded
% (convert point on sphere to point on ellipsoid)
% compute radius of curvature in prime vertical
N_phi = semi_major_axis/sqrt(1-esq*sin_phi_old*sin_phi_old);
P = N_phi*cos_phi_old;                        % P is distance from Z axis
init_z = N_phi*(1-esq)*sin_phi_old;           % initial estimate for z 
init_x = P*cos_lambda_old;                    % initial estimate for x
init_y = P*sin_lambda_old;                    % initial estimate for y 
% End of simplified FRGEOD.M unfolded

% Initialize Old_Sum as the squared residual for the initial receiver 
% position estimate
m = length(pseudoranges);            % find number of SVs
for t = 1:m
    sat_clock(t) = tcorr(t)*vlight;  % compute tcorrs influence on SV clocks
    cal_one_way(1,t) = norm(XS(:,t)-[init_x,init_y,init_z]'); % find distance to initial position
    one_way_res(1,t) = pseudoranges(t,1)-cal_one_way(1,t)+sat_clock(t);
end;
resid_t = one_way_res(1,:)-one_way_res(1,1);    % compute residual for initial position
Old_Sum = resid_t*resid_t';                     % initialize Old_Sum for initial position

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main part of code: search for Least-Squares Sum %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iter = 1:iterations                           % begining of the main loop 
    sin_phi0 = sin_phi_old;                       % update sin_phi0 and cos_phi0
    cos_phi0 = cos_phi_old;                       %  
    for b = 1:ndiv,                               % begining of slicing loop (inclination)
        psi_index=iter*ndiv+b-1;                  % psi_index=iter when ndiv=1
        for alpha_index = 1:alpha_counter,        % begining of slicing loop (azimuth)
            % update grid point using spherical triangle
            sin_phi2 = sin_phi0*cos_psi(psi_index)...
                     +cos_phi0*sin_psi(psi_index)*cos_alpha(alpha_index);
            cos_phi2 = sqrt(1-sin_phi2*sin_phi2); % compute sin and cos for new possible phi
            if cos_phi2 == 0                      % check to avoid divide by zero
                sin_dlambda = 0;                  % if cos_phi2 = 0 => sin_dlambda = 0
            else                                  % else sin_dlambda is computed like this
                sin_dlambda = sin_alpha(alpha_index)*sin_psi(psi_index)/cos_phi2;% a real division
            end;                                  % end of divide by zero check
            
            % To avoid computation of dlambda using asin and thereafter updating lambda2 by adding 
            % lambda_old and dlambda, we instead compute cos_dlambda by sqrt and update sin_lambda2
            % using the fact that sin(A+B)=sin(A)*cos(B)+cos(A)*sin(B). The same is done for cos_lambda2
            cos_dlambda = sqrt(1-sin_dlambda*sin_dlambda);
            sin_lambda2 = sin_lambda_old*cos_dlambda+cos_lambda_old*sin_dlambda;
            cos_lambda2 = cos_lambda_old*cos_dlambda-sin_lambda_old*sin_dlambda;
            
            % The new posible (phi,lambda) coordinates (which are really only known as the sine
            % and cosine of the angles) are converted to Carthesian koordinates (X,Y,Z) using 
            % an simplified unfolded version of FRGEOD.M
            N_phi = semi_major_axis/sqrt(1-esq*sin_phi2*sin_phi2);
            P = N_phi*cos_phi2;
            XR(3,1) = N_phi*(1-esq)* sin_phi2; % XR is a column vector with receiver coordinates
            XR(1,1) = P*cos_lambda2;           % XR(1,1) ~ x, XR(2,1) ~ y, XR(3,1) ~ z
            XR(2,1) = P*sin_lambda2;           % 
            %     End of simplified FRGEOD.M unfolded
            
            for t = 1:m                              % satellite loop
                cal_one_way(1,t) = norm(XS(:,t)-XR); % distance to one of the SVs
                % compute the difference between the pseudorange and the new posible position
                one_way_res(1,t) = pseudoranges(t,1)-cal_one_way(1,t)+sat_clock(t);
            end;                                     % end of satellite loop
            resid_t = one_way_res(1,:)-one_way_res(1,1);% difference between the first 
                                                     % distance and the others
            New_Sum = resid_t*resid_t';              % compute sum of squares of the residuals
            if New_Sum < Old_Sum,                    % test if this sum is the least so far
                Old_Sum = New_Sum;                   % if so, save it
                sin_phi_old = sin_phi2;              % and save (phi,lambda) for next iteration
                cos_phi_old = cos_phi2;              %
                sin_lambda_old = sin_lambda2;        %
                cos_lambda_old = cos_lambda2;        %
            end;                                     % end of test
        end                                          % end of slicing loop (azimuth) 
    end                                              % end of slicing loop (inclination) 
end                                                  % end of main loop (iter)

% The main loop is finished and now the result must be corrected for the 
% rotation of the Earth, while the signals travelled through space:                

corr = (Omegae_dot*mean(pseudoranges(:,1))*inv_vlight)*inv_dtr;  % correction angle in radians
if sign(cos_lambda_old) == 1 & sign(sin_lambda_old) == 1
    lambda_old = asin(sin_lambda_old)*inv_dtr;
end;
if sign(cos_lambda_old) == 1 & sign(sin_lambda_old) == -1
    lambda_old = 360+asin(sin_lambda_old)*inv_dtr;
end;
if sign(cos_lambda_old) == -1 
    lambda_old = 180-asin(sin_lambda_old)*inv_dtr;
end;

lambda = mod(lambda_old-corr,360);  % correct lambda (only 
% lambda is corrected)
phi = asin(sin_phi_old)*inv_dtr;    % compute final phi
%%%%%%%%%%%%%%%%%% end recpos_s.m %%%%%%%%%%%%%%%%%%%%%













