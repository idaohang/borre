function [x,P] = k_update(x,P,A,b,R)
%K_UPDATE   Kalman update, one measurement per call
%	         Observation covariance R

%Written by Kai Borre
%November 24, 1996

omc = b-A*x;
AP = A*P;
innovation_variance = AP*A'+R;
K = AP'/innovation_variance;
x = x+K*omc;
P = P-K*AP;
%%%%%%%% end k_update.m  %%%%%%%%%%%%%%%%%%
