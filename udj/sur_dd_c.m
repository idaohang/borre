%SUR_DD_C  Surface plot of the covariance matrix
%          for double differenced observations

%Copyright (c) by Kai Borre
%$Revision: 1.0$  $Date1999/11/03  $

subplot(2,2,1); surf(dd_cov(2,8)); title('dd\_cov(2,8)')
subplot(2,2,2); surf(dd_cov(3,8)); title('dd\_cov(3,8)')
subplot(2,2,3); surf(dd_cov(4,8)); title('dd\_cov(4,8)')
subplot(2,2,4); surf(dd_cov(5,8)); title('dd\_cov(5,8)')
%%%%%%%%%%%%%%%%%%%%%%%  sur_dd_c  %%%%%%%%%%%%%%%%%%
