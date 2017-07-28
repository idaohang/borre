function [X,P,K,innovation_variance] = k_updatx(X,P,H,Y,R,Q)
%  K_UPDATX   Kalman update, one measurement per call
%	      Allows for system covariance Q
%	      Allows for observation covariance R

%  Written by Kai Borre
%  December 18, 1996

     omc = Y-H'*X;
     P = P + Q;
     HP = H'*P;
     innovation_variance = HP*H+R;
     K = HP'/innovation_variance;
     X = X+K*omc;
     P = P-K*HP;

%%%%%%%% end k_updatx.m  %%%%%%%%%%%%%%%%%%
