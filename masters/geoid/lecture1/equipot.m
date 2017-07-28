%EQUIPOT   It is illustrative to realize that the shape 
%          of the equipotential surfaces of V depend on
%          the dominance of the appropriate terms in the
%          series. An equatorial bulge is caused by the
%          presence of the (2,0) term; the first figure
%          plot teh associated Legenre polynomials P_{2,0}
%          and P_{3,0}. Note that they have 2 and three 
%          zeros, respectively. The modelling of
%          the bulge by the zonal harmonic P_{2,0} is shown
%          in the second figure. The pear-shapedness is 
%          caused by the presence of the (3,0) term, confer
%          the last figure

% Written by Kai Borre
% March 19, 2000

xo = 0:5:180;
x = xo*pi/180;
P2 = (3*(cos(x)).^2-1)/2;
P3 = (5*(cos(x)).^3-3*cos(x))/2;

root2 = fzero('(3*(cos(x)).^2-1)/2',[0 1]);
fprintf('\n First root of P2 is %6.3f degrees',root2*180/pi)
root3 = fzero('(5*(cos(x)).^3-3*cos(x))/2',[0 1]);
fprintf('\n First root of P3 is %6.3f degrees\n',root3*180/pi)

figure(1);
plot(xo,P2)
hold on
plot(xo,P3)
plot([0 180],[0 0],'r')
hold off
title('Associated Legendre Functions {\itP}_{20} and {\itP}_{30}')
text(30,0.7,'{\itP}_{20}')
text(20,0.3,'{\itP}_{30}')

figure(2);
A20 = -0.2;
A30 = 0.15;

add2 = A20*P2;
add2x = add2.*cos(x);
add2y = add2.*sin(x);
add3 = A30*P3;
add3x = add3.*cos(x); 
add3y = add3.*sin(x); 

xx = linspace(0,180,37);
xxo = xx*pi/180;
plot(sin(xxo),cos(xxo),'g', ...
   -sin(fliplr(xxo)),cos(fliplr(xxo)),'g',...
         [-1 1],[0 0],'c',[0 0],[-1 1],'c',...
         [0 sin(root2)],[0 cos(root2)],'c');
text(.03,.2,[num2str(root2*180/pi,3) '{\it\circ}'])
hold on
pl2 = plot(sin(xxo)+add2y,cos(xxo)+add2x,'r',...
                   -sin(fliplr(xxo))-fliplr(add2y),...
                    cos(fliplr(xxo))+fliplr(add2x),'r');
set(pl2,'linewidth',2)
axis equal
hold off
title('Equatorial Bulge')

figure(3);
plot(sin(xxo),cos(xxo),'g', ...
         -sin(fliplr(xxo)),cos(fliplr(xxo)),'g',...
         [-1 1],[0 0],'c',[0 0],[-1 1],'c',...
         [0 sin(root3)],[0 cos(root3)],'c');
text(.03,.3,[num2str(root3*180/pi,3) '{\it\circ}'])
hold on
pl3 = plot(sin(xxo)+add3y,cos(xxo)+add3x,'r',... 
              -sin(fliplr(xxo))-fliplr(add3y),...
                   cos(fliplr(xxo))+fliplr(add3x),'r');
set(pl3,'linewidth',2)
axis equal
hold off
title('Pear-Shapedness')

%%%%%%%%%%%%%%%%%%%%%%%%% end equipot.m  %%%%%%%%%%%%%%%%%%%
