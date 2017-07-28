%L2     Script for computing "En rektangulær grund"

%Written by Kai Borre
%October 6, 1999

format bank
f = [128.1 62.5;
   128.1 62.7;
   128.2 62.6;
   128.0 62.6;
   128.1 62.5];
mf = mean(f); % mf(1)=e, mf(2)=g
d = norm(mf)
a = mf(1)*mf(2)

syms a d e g real;
d = sqrt(e^2+g^2);
a = e*g;
B = jacobian([d; a],[e; g]);
subs(B,{sym('1/(e^2+g^2)^(1/2)')},{sym('1/d')});

syms sigma_e sigma_eg sigma_g;
Sigma = [sigma_e sigma_eg; sigma_eg sigma_g];
Sigma_F = B*Sigma*B';
simplify(Sigma_F);
covf = cov(f);
subs(Sigma_F,{e,g,sigma_e,sigma_eg,sigma_g}, ...
             {mf(1),mf(2),covf(1,1),covf(1,2),covf(2,2)});
sigma_d = sqrt(Sigma_F(1,1));
format short
subs(sigma_d,{e,g,sigma_e,sigma_eg,sigma_g},...
              {mf(1),mf(2),covf(1,1),covf(1,2),covf(2,2)})
sigma_a = sqrt(Sigma_F(2,2));
subs(sigma_a,{e,g,sigma_e,sigma_eg,sigma_g},...
              {mf(1),mf(2),covf(1,1),covf(1,2),covf(2,2)})
format bank
%%%%%%%%%%%%%%%%%% end l2.m  %%%%%%%%%%%%%%%%%%