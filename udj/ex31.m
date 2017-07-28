% Eksempel 3.1 i Mindste kvadraters princip, side 63--65

% Kai Borre, October 13, 1999

given = [270.71 170.71;
         100.00 100.00;
         100.00 241.42];
X_0 = [170.71 170.71];
obs = [100.01; 100.02; 100.03];
deltaY = given(:,1)-X_0(1,1);
deltaX = given(:,2)-X_0(1,2);
alpha = atan2(deltaY,deltaX);
A = [-sin(alpha) -cos(alpha)];
b = obs-sqrt((deltaY(:,1)).^2+(deltaX(:,1)).^2);
x_hat = A\b;
X_0 = X_0+x_hat';
fprintf('\n Y_P: %6.3f  X_P: %6.3f\n',X_0)
r_hat = A*x_hat-b;
sigma2_0 = r_hat'*r_hat/(size(A,1)-size(x_hat,1));
Sigma_xhat = sigma2_0*inv(A'*A);
fprintf('\n sigma Y_P: %3.3f sigma X_P: %3.3f\n\n',... 
   sqrt(Sigma_xhat(1,1)),sqrt(Sigma_xhat(2,2)))

% New point of development
fprintf('\n With new point of development')
i = 0;
test = 10;
X_0 = [100  100]; % far from the solution
while test > 1.e-3
   i = i+1;
   fprintf('\n Iteration %3.0f\n',i)
   deltaY = given(:,1)-X_0(1,1);
   deltaX = given(:,2)-X_0(1,2);
   alpha = atan2(deltaY,deltaX);
   A = [-sin(alpha) -cos(alpha)];
   b = obs-sqrt((deltaY(:,1)).^2+(deltaX(:,1)).^2);
   x_hat = A\b;
   test = norm(x_hat);
   X_0 = X_0+x_hat';
   fprintf(' Y_P: %6.3f  X_P: %6.3f\n',X_0)
   r_hat = A*x_hat-b;
   sigma2_0 = r_hat'*r_hat/(size(A,1)-size(x_hat,1));
   Sigma_xhat = sigma2_0*inv(A'*A);
   fprintf('\n sigma Y_P: %3.3f sigma X_P: %3.3f\n',... 
      sqrt(Sigma_xhat(1,1)),sqrt(Sigma_xhat(2,2)))
end

%%%%%%%%%%%%%%%%%%%%% end ex31.m  %%%%%%%%%%%%%%%%%

