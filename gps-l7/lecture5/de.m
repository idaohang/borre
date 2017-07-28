%DE     Analysis of the filter matrix of the four
%       observation filter N1 and N2. Transformation to
%       the wide lane ambiguity Nw and N1

%Kai Borre, July 22, 1995
%Copyright (c) by Kai Borre
%$Revision: 1.0$  $Date:1998/11/08  $

e = exist('de1.ps');
if e ~= 0, delete de1.ps, end
e = exist('de2.ps');     
if e ~= 0, delete de2.ps, end

%  Some constants
c0 = 299792458;	% velocity of light in m/s
f1 = 154*10.23E6; % L1 frequency in Hz
f2 = 120*10.23E6; % L2 frequency in Hz
lambda1 = c0/f1;	% .19029367  m
lambda2 = c0/f2;	% .244210213 m
%  Initialization of A matrix in filter
A = [ones(4,2) zeros(4,2)];
A(1,2) = 1;
A(3,2) = (f1/f2)^2;
A(2,2) = -1;
A(4,2) = -(f1/f2)^2;
A(2,3) = lambda1;
A(4,4) = lambda2;
A
sigma_1Phi = 0.002;   % std. dev. of phase obsv.
k = 154;		          % sigma_i,P - sigma_i,Phi
v = [ k^2; 1; (k*f1/f2)^2; (f1/f2)^2];
Sigma_b = sigma_1Phi^2*diag(v);
Sigma_x = inv(A)*Sigma_b*(inv(A))'
for i = 1:4
   for j = 1:4
      corr_x(i,j) = Sigma_x(i,j)/sqrt(Sigma_x(i,i)*Sigma_x(j,j));
   end
end
corr_x

[V,D] = eig(Sigma_x(3:4,3:4));
[lambda,m] = sort(diag(D));
V = V(:,m);
if any(any(V)) == 1
   alpha = atan2(V(1,2),V(2,2)); 
end
rot = [cos(alpha) -sin(alpha);sin(alpha) cos(alpha)];
alpha = alpha*180/pi	       % convertion to degrees
t = linspace(0,2*pi,50);
a = sqrt(lambda(2))
b = sqrt(lambda(1))
pl = [a*sin(t);b*cos(t)];
for t = 1:50
   current = rot*pl(:,t);
   curve(1:2,t) = current; 
end

fig1 = figure;
axis([-12 12 -12 12])
%      set(gca,'XLim',[-10 10]);
plot(curve(1,1:50),curve(2,1:50))
xlabel('\itN\rm_2','FontSize',20)
ylabel('\itN\rm_1','FontSize',20);
set(gca,'FontSize',20);
grid  on
print de1 -deps

%% transformation of observations to yield Nw and N1
T = eye(4);
T(4,4) = 0;
T(4,3) = 1;
T(3,4) = -1;

Sigma_z = T*inv(A)*Sigma_b*(inv(A))'*T'
for i = 1:4
   for j = 1:4
      corr_z(i,j) = Sigma_z(i,j)/sqrt(Sigma_z(i,i)*Sigma_z(j,j));
   end
end
corr_z

[V,D] = eig(Sigma_z(3:4,3:4));
[lambda,m] = sort(diag(D));
V = V(:,m);
if any(any(V)) == 1
   alpha = atan2(V(2,2),V(1,2)); 
end
rot = [cos(alpha) -sin(alpha);sin(alpha) cos(alpha)];
alpha = alpha*180/pi;		  % convertion to degrees
if alpha < 0
   alpha = alpha+180; 
end;
alpha
t = linspace(0,2*pi,50);
a = sqrt(lambda(2))
b = sqrt(lambda(1))
pl = [a*sin(t);b*cos(t)];
for t = 1:50
   current = rot*pl(:,t); 
   curve(1:2,t) = current;
end

fig2 = figure;
plot(curve(1,1:50),curve(2,1:50))
axis([-12 12 -12 12]);
xlabel('\itN\rm_1-\itN\rm_2','FontSize',20);
ylabel('\itN\rm_1','FontSize',20)
set(gca,'FontSize',20);
grid
print de2 -deps
%%%% end de.m %%%%%%%%%%%%%%%%%%%%%
