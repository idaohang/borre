function svdshow(A)
%SVDSHOW Graphical demonstration of matrix singular values.
%   SVDSHOW(A), where A is a 2-by-2 matrix, presents a graphical
%   experiment involving the effect of the mapping induced by A on
%   the unit circle.  Use the mouse to move the vector x and follow
%   the resulting A*x and A*y.  When A*x is perpendicular to A*y,
%   then x and y are right singular vectors of A and the length of
%   A*x and A*y are the corresponding singular values.
%   SVDSHOW, by default, uses
%      A = [1/4 3/4; 1 1/2];
%   Other interesting examples are:
%      A = eye(2,2);
%      A = [0 1; 1 0];
%      A = [0 1; -1 0];
%      A = [1/2 1; -1/4 -1/2]
%      A = [3/2 1; -1/4 1/2]
%      A = randn(2,2);
%
%   See also EIGSHOW.

%   Copyright (c) 1993-97 by The MathWorks, Inc.
%   $Revision: $  $Date: $

if nargin == 0;
   A = [1/4 3/4; 1 1/2];
   svdshowinit(A);
elseif all(size(A) == 2)
   svdshowinit(A);
elseif isequal(A,'action')
   svdshowaction
else
   error('Matrix must be 2-by-2.')
end

%------------------

function svdshowinit(A)
clf reset
s = 1.1*max(1,norm(A));
axis([-s s -s s])
axis square
title(sprintf('A = [ %6.3f  %6.3f; %6.3f  %6.3f ]',A'))
xlabel('Make A*x perpendicular to A*y')
xcolor = [0 .6 0];
Axcolor = [0 0 .8];
h.x = eigshowvec([1 0]','x',xcolor);
h.Ax = eigshowvec(A(:,1),'Ax',Axcolor);
h.y = eigshowvec([0 1]','y',xcolor);
h.Ay = eigshowvec(A(:,2),'Ay',Axcolor);
set(gcf,'userdata',A);
set(gca,'userdata',h);
set(gcf,'windowbuttondownfcn', ...
   'svdshow(''action'');set(gcf,''windowbuttonmotionfcn'',''svdshow(''''action'''')'')')
set(gcf,'windowbuttonupfcn', ...
  'set(gcf,''windowbuttonmotionfcn'','''')')

%------------------

function svdshowaction
A = get(gcf,'userdata');
h = get(gca,'userdata');
pt = get(gca,'currentpoint');
x = pt(1,1:2)';
x = x/norm(x);
y = [-x(2) x(1)]';
eigshowact(x,h.x);
eigshowact(y,h.y);
eigshowact(A*x,h.Ax);
eigshowact(A*y,h.Ay);

%------------------

function h = eigshowvec(x,t,color)
h.mark = line(x(1),x(2),'marker','.','erase','none','color',color);
h.line = line([0 x(1)],[0 x(2)],'erase','xor','color',color);
h.text = text(x(1)/2,x(2)/2,t,'fontsize',12,'erase','xor','color',color);


%------------------

function eigshowact(x,h)
set(h.mark,'xdata',x(1),'ydata',x(2));
set(h.line,'xdata',[0 x(1)],'ydata',[0 x(2)]);
set(h.text,'pos',x/2);
