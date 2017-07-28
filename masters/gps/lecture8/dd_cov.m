function M = dd_cov(r,s);
%DD_COV  Computes the covariance matrix for double 
%        differenced observations between r receivers
%        and s satellites

%Written by Kai Borre
%January 15, 1998

r = r-1;
s = s-1;
R = [ones(r,1) -eye(r)];
S = [-ones(s,1) eye(s)];
Dd = kron(R,S);
Sigma_d = Dd*Dd';
M = (r+1)*(s+1)*inv(Sigma_d);
%%%%%%%%%%%% end dd_cov  %%%%%%%%%%%%%%%%%%%%%%
