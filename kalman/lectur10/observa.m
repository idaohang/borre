%OBSERVA  Test for observability, and filtering of double
%         differenced pseudoranges for static receivers
%         on L1 and L2.

%Written by Kai Borre and Christian C. Tiberius, December 19, 1998
%Copyright (c) by Kai Borre and Christian C. Tiberius
%$Revision: 1.0$  $Date:1998/12/20  $

% units in meters and seconds
f1 = 154*10.23e6;
f2 = 120*10.23e6;
beta = (f1/f2)^2;
dt = 20;

% Setup of model
% For DD data a first order model for rho works well.

% Take rho_ddot and I_ddot as noise terms, model them as a white
% noise random process with zero mean

% Coefficient matrix for dual frequency pseudorange data
A = [1 0 1    0;
     1 0 beta 0];

% State vector x = (rho, rho_dot, I, I_dot)  

% State transition matrix
% rho: first order model, rho_ddot = noise
% iono: first order model, I_ddot = noise
F = [ 1 dt  0  0;
      0  1  0  0;
      0  0  1 dt;
      0  0  0  1];

% Observation covariance matrix R.
R = 0.1*eye(2);

% State covariance matrix Q.
% The spectral density of rho_ddot [m2/s3]
s_r = 0.36e-8; 
% The spectral density of I_ddot [m2/s3] 
s_I = 4.e-8; 
% Covariance matrix for first order model
Q1 = [1/3*dt^3 1/2*dt^2;
      1/2*dt^2      dt];
Q = [Q1*s_r     zeros(2,2);
     zeros(2,2)    Q1*s_I];

data = [...
   -3637.5195    -3638.0597;
   -3632.7182    -3632.7425;
   -3626.8818    -3626.9856;
   -3621.1026    -3621.2142;
   -3615.7202    -3615.7936;
   -3610.5369    -3610.0575;
   -3604.4218    -3604.5555;
   -3598.4928    -3598.8831;
   -3593.2302    -3593.3801;
   -3588.0704    -3588.0644;
   -3582.1526    -3582.3553;
   -3576.5273    -3576.9803;
   -3571.0656    -3571.4433;
   -3565.5150    -3566.1544;
   -3560.4823    -3560.6687;
   -3554.9714    -3555.1840;
   -3549.4712    -3549.3614;
   -3544.1538    -3544.1706;
   -3538.7740    -3538.6484;
   -3533.6450    -3533.3713;
   -3528.0363    -3528.0795;
   -3522.5657    -3522.7747;
   -3517.2219    -3517.6254;
   -3511.9404    -3512.1836;
   -3506.7126    -3506.9548;
   -3501.4295    -3501.7704;
   -3496.2624    -3496.4960;
   -3490.9270    -3491.0454;
   -3485.7241    -3485.9378;
   -3480.6678    -3480.8522;
   -3475.4775    -3475.7545;
   -3470.2300    -3470.5145;
   -3465.1446    -3465.3739;
   -3460.0277    -3460.1023;
   -3455.3969    -3454.8953;
   -3450.3986    -3449.7804;
   -3445.0875    -3445.1052;
   -3440.2709    -3440.0112;
   -3435.2351    -3435.1340;
   -3430.3882    -3430.3117;
   -3425.3985    -3425.4593;
   -3420.2264    -3420.4624;
   -3415.5292    -3415.4396;
   -3410.5864    -3410.4398;
   -3405.5588    -3405.6221;
   -3400.6556    -3400.9731;
   -3395.8556    -3396.0450;
   -3391.1911    -3391.2822;
   -3386.4275    -3386.4587;
   -3381.6682    -3381.7363;
   -3377.1247    -3376.9772;
   -3372.3359    -3372.3419;
   -3367.7623    -3368.0887;
   -3363.1993    -3363.4835;
   -3358.7424    -3358.5383;
   -3354.2484    -3353.9933;
   -3349.5812    -3349.3953;
   -3344.7886    -3345.0389;
   -3340.1925    -3340.4284;
   -3336.0100    -3335.9099;
   -3331.5214    -3331.5163;
   -3326.9547    -3327.0110;
   -3322.1351    -3322.5201;
   -3317.5005    -3317.9568;
   -3313.2910    -3313.5827;
   -3309.0862    -3308.8761;
   -3304.7269    -3304.6205;
   -3300.4191    -3300.1367;
   -3295.9053    -3296.0240;
   -3291.3235    -3291.7675;
   -3287.1127    -3287.4877;
   -3283.3221    -3283.3899;
   -3279.2166    -3279.2669;
   -3275.4575    -3275.2462;
   -3271.4580    -3271.2203;
   -3267.1856    -3266.8562;
   -3263.1203    -3262.8688;
   -3259.1539    -3258.7848;
   -3254.8985    -3254.7603;
   -3250.9048    -3250.6716;
   -3246.9010    -3246.8120;
   -3242.6589    -3242.8115;
   -3238.6958    -3238.9359;
   -3234.7311    -3234.8076;
   -3230.7014    -3230.9406;
   -3226.9668    -3226.7495;
   -3223.1473    -3222.8679;
   -3218.9925    -3219.3783;
   -3215.2833    -3215.4949;
   -3211.6432    -3211.8174];

% Initialization:
% The first estimate is based on 2 epochs. The coefficient matrix
% is Ai, covariance matrix Qi, and data-vector yi 
yi = [data(1,:)'; zeros(4,1); data(2,:)'];
Ai = [ A          zeros(2,4);
      -F          eye(4)    ;
      zeros(2,4)  A        ];
deficiency = size(Ai,1)-rank(Ai);
if deficiency == 0
   fprintf('\nThe system is observable after 2 epochs\n') 
end

Qi = [R          zeros(2,4) zeros(2,2);
      zeros(4,2) Q          zeros(4,2);
      zeros(2,2) zeros(2,4) R        ];

% Compute least-squares solution, 
Pi = inv(Ai'*inv(Qi)*Ai);
xi = Pi*Ai'*inv(Qi)*yi;

% Extract last state from this solution
x = xi(5:8);
P = Pi(5:8,5:8);
x_acc = [];
P_acc = [];

for i = 2:size(data,1)
   [x,P,K,innovar] = k_updatf(x,P,A,data(i,:)',R,Q,F);
   x_acc = [x_acc x];
   P_acc = [P_acc P];
end

t = 2:i;
figure(1);
plot(t,x_acc(1,:))
title('Range, [m]')

figure(2);
plot(t,x_acc(2,:))
title('Range Rate, [m/s]')

figure(3);
plot(t,x_acc(3,:))
title('Ionospheric Delay, [m]')

figure(4);
plot(t,x_acc(4,:))
title('Ionospheric Delay Rate, [m/s]')

figure(5);
plot(P_acc(1,1:4:4*(i-2)))
title('Standard Deviation of Range, [m]')
%%%%%%%%%%%%%%%%% end observa.m %%%%%%
