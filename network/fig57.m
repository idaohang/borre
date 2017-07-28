%FIG57   Plot of Weyl's asymptotic theorem for eigenvalues

%Kai Borre April 18, 1998
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

x = 0:1000;
eigs = x/(4*pi)+sqrt(x)/pi;
% exact counting of the integers
eig2 = zeros(1,1001);
for n = 1:335
   for m = 1:350
      eig = (m^2+n^2)*pi^2;
      ind = round(eig);
      if ind < 1000
         eig2(1,ind) = eig2(1,ind)+1;
      end
   end
end
eigacc = cumsum(eig2);
eigacc = eigacc+2; %the two zero eigenvalues!

plot(x,eigs)
hold on
plot(x,eigacc)
hold off
ylabel('\itN\rm(\lambda)','VerticalAlignment','Bottom','FontSize',18)
xlabel('\lambda','FontSize',18)
set(gca,'FontSize',18)
print -deps fig57
%%%%%%%%%%%%%%%%% end fig57.m %%%%%%%%%%%%%