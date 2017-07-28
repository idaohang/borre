function ex1611(sb,sd)
%EX1611   Discrete random ramp for estimating receiver 
%              clock offset and drift.
%	           Covariance for system noise Q,
%              and observation noise R.
%              The observation series b represents receiver
%              clock offsets (ns) from a steered clock.

%Kai Borre 12-05-98
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 1998/12/05  $

if nargin == 0
    sb =  1.e-11;   % sigma of short term stability 
    sd =  1.e-10;   % sigma of long term stability 
end
b = [...
      -406.05 -331.08 -266.85	 -210.83  -170.77 ...
      -134.14  -96.63  -65.40	  -40.09   -23.08 ...
       -17.42  -17.54  157.41	  331.35   493.32 ...
       642.38  752.80  848.02	  925.23   974.17 ...
       980.48  972.66  929.62	  854.92   734.40 ...
       626.17  541.18  475.70	  421.21   374.44 ...
       322.89  254.76  226.38	  190.00   157.51 ...
       139.05  106.41   75.90	   48.49    17.22 ...
       -19.32  -48.17  -78.52	 -113.18  -145.47 ...
      -170.82 -190.06 -207.77	 -218.61  -226.48 ...
      -237.92 -233.96 -223.53	 -210.70  -196.02 ...
      -186.03 -224.64 -246.42	 -238.71  -193.83 ...
      -172.58 -161.17 -152.88	 -153.72  -156.46 ...
      -164.55 -174.27 -173.64	 -158.99  -152.54 ...
      -145.98 -135.08 -120.72	 -119.30  -112.49 ...
       -88.95  -64.95  -49.72	  -25.61    -0.25];

N = size(b,2);
dt = 20*1.e9;                     %epoch interval in nanoseconds
F = [1 dt;0 1];
A = [1 dt];
R = 10;                           %observation error variance
Q = [sb^2+(sd*dt)^2 sd^2*dt;
        sd^2*dt     sd^2  ];      %system error covariance

% Initial conditions
x_minus = [b(1); 1.e-8];
P_minus = [1.e4 0;0 1.e-16]; 
x_acc = [];

%Kalman filtering
for i = 0:N-1
   [x_plus, P_plus, K, innovation_var] = ...
                 k_updatx(x_minus, P_minus, A, b(i+1), R, Q);
   x_minus = x_plus;
   P_minus = F*P_plus*F'+Q;
   x_acc = [x_acc x_minus];
end

t = 0:N-1;
comp = A*x_acc;

figure(1);
plot(t,b(1:N),'mo',t,comp(1:N),'g:');
title('Receiver clock offset','Fontsize',12)
ylabel('[ns]')
set(gca,'Fontsize',12);
print ex16111 -deps

figure(2);
plot(x_acc(1,:)); 
title('Filtered clock offset')
ylabel('[ns]','Fontsize',12)
set(gca,'Fontsize',12);
print ex16112 -deps

figure(3);
plot(x_acc(2,:)); 
title('Filtered clock drift')
ylabel('[ns]','Fontsize',12)
set(gca,'Fontsize',12);
print ex16113 -deps

figure(4);
plot(b-comp); 
title('Observed minus filtered clock offset')
ylabel('Residual [ns]','Fontsize',12)
set(gca,'Fontsize',12);
print ex16114 -deps
%%%%%%%% end ex1611.m  %%%%%%%%%%%%%%%%%%
