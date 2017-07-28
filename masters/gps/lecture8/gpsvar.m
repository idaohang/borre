%GPSVAR Estimation of standard deviations for N1, N2,
%       and N1-N2 when using a least-squares procedure
%       in batch mode. There are 2 code and 2 phase
%       observations, and I is put equal to zero.

%Kai Borre, March 26, 1996
%Copyright (c) by Kai Borre
%$Revision: 1.0$  $Date:1998/11/08  $

% Initial computations of constants
c = 299792458;     % vacuum speed of light in m/s
f1 = 154*10.23E6;  % L1 frequency in Hz
f2 = 120*10.23E6;  % L2 frequency in Hz
lambda1 = c/f1;    % wavelength on L1:  .19029367  in m
lambda2 = c/f2;    % wavelength on L2:  .244210213 in m

e = exist('gpsvar.eps');
if e ~= 0, delete gpsvar.eps, end;

% The model looks like
%  b = [P1; Phi1; P2; Phi2];
%  x = [rho; N1; N2];
% Definition of A matrix
A = [ones(4,1) zeros(4,2)];
A(2,2) = lambda1;
A(4,3) = lambda2;
% Transformation matrix for wide lane
Z1 = eye(2);
Z1(2,1) = -1;
sd = [0.3 0.005 0.3 0.005];		   % standard deviations
W = diag([1./sd.^2]);		         % weight matrix
Wa = diag([1./sd.^2 -1/trace(W)]);	% augmented weight matrix
ev = ones(4,1);		            	% sum does not function
A1 = [A; ev'*W*A];  % A augmented with a fictituous obsv. eq.
A_aug = [];
sigmaN1 = [];
sigmaN2 = [];
sigmaN1_N2 = [];

for epoch = 1:100
   A_aug = [A_aug; A1];
   N = A_aug'*kron(eye(epoch),Wa)*A_aug;
   S = inv(N(2:3,2:3));
   sigmaN1 = [sigmaN1 sqrt(S(1,1))];
   sigmaN2 = [sigmaN2 sqrt(S(2,2))];
   ST = Z1*S*Z1';         % transformation for wide lane
   sigmaN1_N2 = [sigmaN1_N2 sqrt(ST(2,2))];
end

t = 1:100;
figure(1);
hold on
pl1 = plot(t,sigmaN1,'o',t,sigmaN2,'x',t,sigmaN1_N2,'*');
ylabel('Standard deviation [m]','FontSize',16)
xlabel('Number of epochs','FontSize',16)
set(pl1,'Markersize',3);
hold off
set(gca,'Fontsize',16)
text(20,0.32,'\itN\rm_1','FontSize',16)
text(7,0.23,'\itN\rm_2','FontSize',16)
text(12,0.13,'\itN\rm_1-\itN\rm_2','FontSize',16)
print gpsvar -deps
%%%%%%%%%% end gpsvar.m %%%%%%%%%%%%


