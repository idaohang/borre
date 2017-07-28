function [x,P,K,innovation_variance] = k_updatf(x,P,A,b,R,Q,F)
%K_UPDATF   Kalman update, one measurement per call
%  	      Allows for system covariance Q
%	         Allows for observation covariance R

%Written by Kai Borre
%January 27, 1998

x = F*x;
omc = b-A*x;
P = F*P*F'+Q;
AP = A*P;
innovation_variance = AP*A'+R;
K = AP'/innovation_variance;
x = x+K*omc;
P = P-K*AP;
%%%%%%%% end k_updatf.m  %%%%%%%%%%%%%%%%%%
