function [x,P,K,innovation_variance] = k_updatx(x,P,A,b,R,Q)
%K_UPDATX   Kalman update, one observation per call.
%  	                 System covariance Q, and observation covariance R

%Kai Borre, December 18, 1996
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date:1998/12/08  $

omc = b-A*x;
P = P + Q;
AP = A*P;
innovation_variance = AP*A'+R;
K = AP'/innovation_variance;
x = x+K*omc;
P = P-K*AP;
%%%%%%%% end k_updatx.m  %%%%%%%%%%%%%%%%%%
