% DATA_ANA  Error analysis of GPS observables by making linear 
%            combinations. We use data from one receiver to one 
%            satellite.
%            Idea by C. Tiberius

% Written by Kai Borre, February, 21, 2002
% Copyright by Kai Borre 

% Observations from one satellite in a RINEX file are read into a matrix.
% The extension gives the PRN. Each row corresponds to the observations 
% from one epoch; the columns contain the follwing information
% 1         2      3       4     5      6 
% C/A=P1    Phi1   Phi2    P2    S1     S2
% [m]       [cyc]  [cyc]   [m]   [ ]    [ ] 
% S1 and S2 are not used; column 5 alternatively may contain P1

load trm1.001;
x = trm1;
[N hlp1] = size(x);            % determine sample size N

% convert phase to metres; all observations then in [m]
v_light = 299792458;         % velocity of light [m/s]
f1 = 154*10.23e6;            % frequency of L1
f2 = 120*10.23e6;            % frequency of L2
lambda1 = v_light/f1;        % L1 wavelength [m]
lambda2 = v_light/f2;        % L2 wavelength [m]
alpha = (f1/f2)^2;           % alpha = (f1/f2)^2 

x(:,2) = x(:,2)*lambda1;
x(:,3) = x(:,3)*lambda2;

% build Mc combination: 1.00*P1 - 4.09*Phi1 + 3.09*Phi2
%   is a function of only the phase ambiguities and should thus be constant

Mc = x(:,1) + (1+alpha)/(1-alpha)*x(:,2) - 2/(1-alpha)*x(:,3);

plot(Mc-mean(Mc));
title('Mc - PRN04');

% build M2 combination: 1.00*C1 - 5.09*L1 + 4.09*L2
%   is a function of only the phase ambiguities and should thus be constant

M2 = x(:,4) + 2*alpha/(1-alpha)*x(:,2) - (1+alpha)/(1-alpha)*x(:,3);


% build Ic combination: 1.55*(P2-C1)
%   is a function of only the ionospheric delay (from code only)

Ic = 1/(alpha-1)*(x(:,4)-x(:,1));

% build Il combination: -1.55*(L2-L1)
%   is a function of only the ionospheric delay and the constant
%   ambiguities (from phase only)

Il = -1/(alpha-1)*(x(:,3)-x(:,2));

%%%%%%%%%%%%%%%%%%% end err_src.m  %%%%%%%%%%%%%%%%%%%%%%%%%