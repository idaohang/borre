% Solution to Eksempel 6.8
% October 31, 1999

format rat
A = [-1 1 0 0 0 0;
     0 -1 0 0 1 0;
     0 0 0 1 -1 0;
     1 0 0 -1 0 0;
     0 0 0 0 -1 1;
     0 0 1 0 0 -1;
     0 1 -1 0 0 0]
PinvA = pinv(A)
Sigma_xplus = PinvA*PinvA'
f = [-1 0 0 0 0 1]
sigma16 = f*pinv(A'*A)*f'
A_fixed = A(:,[2:5])
N_fixed = A_fixed'*A_fixed
Q_fixed = inv(N_fixed)
g = [-1 0 0 1];
sigma25 = g*Q_fixed*g'
h = [0 -1 0 0 1 0];
sigma25_fixed = h*pinv(A'*A)*h'
fprintf('\n Trace of fixed network:  %8.4f',trace(Q_fixed))
fprintf('\n Trace of free network: %8.4f\n',...
  Sigma_xplus(2,2)+Sigma_xplus(3,3)+Sigma_xplus(4,4)+Sigma_xplus(5,5))
%%%%%%%%%%%%%%%%%%%%%  eks68.m  %%%%%%%%%%%%%%%
