%LEG_POL  Plot of the first 20 Legendre polynomials

%Written by Kai Borre
%March 25, 2000

for n = 1:20 % n is the degree, the order is always equal to m = 0
   tt = 0;
   for x = -1:.01:1
      tt = tt+1;
      [P,Pf] = legen_nm(n-1,0,x);
      L(n,tt) = Pf;
   end
end
xx = linspace(-1,1,201)
plot(xx,L')
title('Legendre Polynomials')

%%%%%%%%%%%%%%%%%%% end leg_pol.m %%%%%%%%%%%%%
