% Example 3.1 in Mindste kvadraters princip, page 63--65
% We define the network as a free network. Next we
% augment the normals, and finally we fix it!

%Copyright (c) Kai Borre
%$Revision: 1.0 $ $Date: 1999/11/06  $

prelim = [270.71 170.71;  % Y,X for 001
          100.00 100.00;  %         002
          100.00 241.42]; %         003
X_0 = [170.71 170.71];    % Y,X for P
       
% We define a least squares problem with 8 unknown coordinates.
% We introduce a G matrix with 3 columns, and include 3 observations.
% So we need at least 2 more observations that the network be 
% defined geometrically. For symmetry reasons we choose 3 additional
% observations (3+3+3 > 8).
obs = [100.01; 100.02; 100.03; 184.785; 141.44; 184.805];
for i = 1:3
   comp(i,1) = norm(prelim(i,:)-X_0);
end
comp(4,1) = norm(prelim(1,:)-prelim(2,:));
comp(5,1) = norm(prelim(2,:)-prelim(3,:));
comp(6,1) = norm(prelim(3,:)-prelim(1,:));
b = obs-comp;  %omc
% The sequence of unknowns: x = (Y1,X1,Y2,X2,Y3,X3,YP,XP)
% We compute the non-zero entries of A. Note that we only
% need two numbers for each observation
for i = 1:3
   s(i) = (prelim(i,1)-X_0(1,1))/comp(i); %sin alpha_1P
   c(i) = (prelim(i,2)-X_0(1,2))/comp(i); %cos alpha_1P
end
s(4) = (prelim(1,1)-prelim(2,1))/comp(4);
s(5) = (prelim(2,1)-prelim(3,1))/comp(5);
s(6) = (prelim(3,1)-prelim(1,1))/comp(6);
c(4) = (prelim(1,2)-prelim(2,2))/comp(4);
c(5) = (prelim(2,2)-prelim(3,2))/comp(5);
c(6) = (prelim(3,2)-prelim(1,2))/comp(6);
A = [ s(1) c(1)    0     0     0     0  -s(1) -c(1);
         0    0  s(2)  c(2)    0     0  -s(2) -c(2);
         0    0    0     0   s(3)  c(3) -s(3) -c(3);
      s(4) c(4) -s(4) -c(4)    0     0      0     0;
         0    0  s(5)  c(5) -s(5) -c(5)     0     0;
      -s(6) -c(6)  0      0  s(6)  c(6)     0     0];     
A0 = A;      
x_plus = pinv(A)*b;
X_plus = X_0+x_plus(7:8,1)'; 
fprintf('\n Y_P_free: %6.4f  X_P_free: %6.4f\n',X_plus)
% One way of removing the singularity of A is by augmenting
% the normals by all columns spanning A's null space
coord0 = [prelim; X_0];
coord = fliplr(coord0);
kk = size(coord,1);
plusminus = eye(2*kk);
for i = 1:2*kk
   if  mod(i,2) ~= 0
      plusminus(i,i) = -1;
   end
end
G = [kron(ones(kk,1),[1;0]), kron(ones(kk,1),[0;1]), ...
       plusminus*reshape(coord',2*kk,1), reshape(coord0',2*kk,1)];         
N_aug = [A'*A G; G' zeros(size(G,2),size(G,2))];
b_aug = [A'*b; zeros(size(G,2),1)];
x_aug = N_aug\b_aug ;
X_aug = X_0+x_aug(7:8,1)';
fprintf('\n Y_P_aug: %6.4f  X_P_aug: %6.4f\n',X_aug)
% We fix the coordinates of points 1, 2, and 3. This is
% achieved by deleting the first 6 columns of A:
A(:,1:6) = [];
% Now A is regular and the solution is obtained by
% using the backslash operator
x_fixed = A\b;
X_fixed = X_0+x_fixed'; 
fprintf('\n Y_P_fixed: %6.4f  X_P_fixed: %6.4f\n',X_fixed)

r_hat = A*x_fixed-b;
sigma2_0 = (norm(r_hat))^2/(size(A,1)-size(x_fixed,1));
Sigma_xhat = sigma2_0*inv(A'*A);
fprintf('\n sigma Y_P: %3.4f sigma X_P: %3.4f\n\n',... 
              sqrt(Sigma_xhat(1,1)),sqrt(Sigma_xhat(2,2)))
% We use an S-transformation as described in
% Landmaaling, page 138 to transform Sigma_xhat 
% into a covariance matrix fixed at points 1 and 2
sg = size(G,2);
S = [zeros(sg,2*kk); 
     -G(sg+1:2*kk,:)*inv(G(1:sg,:)) eye(sg)];
Cov = (norm(b-A0*x_plus))^2*pinv(A0'*A0);
Sigma_xp = S*Cov*S'%;
%%%%%%%%%%%%%%%%%%%%% end ex31free.m  %%%%%%%%%%%%%%%%%

