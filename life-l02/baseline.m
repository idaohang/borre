function bas = baseline(master,obs,sats,time,Eph)
% BASELINE Computation of baseline between master and rover 
%          from pseudoranges alone, and using ordinary 
%          least-squares estimation

%Kai Borre 31-10-2001
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2001/10/31  $

v_light = 299792458;
dtr = pi/180;
m = size(obs,1);  % number of svs
% identify ephemerides columns in Eph
for t = 1:m
    col_Eph(t) = find_eph(Eph,sats(t),time);
end

% preliminary guess for receiver position and receiver clock offset
bas = zeros(4,1);
x = zeros(4,1);
pos = master+bas;
no_iterations = 3; 
for iter = 1:no_iterations
    A = [];
    omc = []; % observed minus computed observation
    for i = 1:m
        k = col_Eph(i);
        tx_RAW = time - obs(i)/v_light;
        t0c = Eph(21,k);
        dt = check_t(tx_RAW-t0c);
        tcorr = (Eph(2,k)*dt + Eph(20,k))*dt + Eph(19,k);
        tx_GPS = tx_RAW-tcorr;
        dt = check_t(tx_GPS-t0c);
        tcorr = (Eph(2,k)*dt + Eph(20,k))*dt + Eph(19,k);
        tx_GPS = tx_RAW-tcorr;
        X = satpos(tx_GPS, Eph(:,k));
        if iter == 1
            traveltime = 0.072;
            Rot_X = X;
            trop = 0;
            rho2 = (norm(X(1:3)-master(1:3),'fro'))^2;
        else
            rho2 = (X(1)-pos(1))^2+(X(2)-pos(2))^2+(X(3)-pos(3))^2;
            traveltime = sqrt(rho2)/v_light;
            Rot_X = e_r_corr(traveltime,X);
            rho2 = (Rot_X(1)-pos(1))^2+(Rot_X(2)-pos(2))^2+(Rot_X(3)-pos(3))^2;
        end
        a =  [-(Rot_X(1)-pos(1))/sqrt(rho2)...
              -(Rot_X(2)-pos(2))/sqrt(rho2) ...
              -(Rot_X(3)-pos(3))/sqrt(rho2) 1];
        % subtraction of pos(4) corrects for receiver clock offset 
        omc = [omc; obs(i)-a*x]; 
        A = [A;a];     
    end % i
    x = A\omc;
    bas = bas-x;
end % iter

%%%%%%%%%%%%%%%%%%%%%  baseline.m  %%%%%%%%%%%%%%%%%%%%%
