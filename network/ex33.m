%EX33   Plot of covariance function for difference of 
%       abscissae between two points r apart and with
%       direction angle phi beetween them.

%Kai Borre May 24, 1999
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2000/12/16 $

h = -2:.05:2;
lh = length(h);
k = -2:.05:2;
lk = length(k);
A = zeros(lh,lk);
B = A;

for r = 1:lh
   for s = 1:lk
      z = h(r)+i*k(s);
      if h(r) ~= 0 | k(s) ~= 0
        A(r,s) = log(abs(z))-0.25*cos(2*angle(z));
        B(r,s) = log(abs(z))+0.25*cos(2*angle(z));
      else
         A(r,s) = NaN;
         B(r,s) = NaN;
      end
   end
end

subplot(2,1,1), meshc(A) %surfc(A);
set(gca,'XTick',[1:(lh-1)/4:lh],'XTickLabel', ...
                    ['-2';'-1';' 0';' 1';' 2'],...
        'YTick',[1:(lk-1)/4:lk],'YTickLabel', ...
                    ['-2';'-1';' 0';' 1';' 2'])
title(['Variance of abscissae difference'... 
      ' between origin and any point (\itX,\itY)'])                
xlabel('\itX')
ylabel('\itY')
           
subplot(2,1,2), surfc(B);
set(gca,'XTick',[1:(lh-1)/4:lh],'XTickLabel', ...
                    ['-2';'-1';' 0';' 1';' 2'],...
        'YTick',[1:(lk-1)/4:lk],'YTickLabel', ...
                    ['-2';'-1';' 0';' 1';' 2'])
title(['Variance of ordinate difference'... 
       ' between origin and any point (\itX,\itY)'])                
xlabel('\itX')
ylabel('\itY')

print -deps ex33
%%%%%%%%%%%%%%%%%%%% end ex33.m  %%%%%%%%%%%%%
