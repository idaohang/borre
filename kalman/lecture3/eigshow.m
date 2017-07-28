function eigshow(A)
%EIGSHOW Graphical demonstration of matrix eigenvalues.
%   EIGSHOW(A), where A is a 2-by-2 matrix, presents a graphical
%   experiment involving the effect of the mapping induced by A on
%   the unit circle.  Use the mouse to move the vector x and follow
%   the resulting A*x.  When A*x is parallel to x, then x is an
%   eigenvector of A and the length of A*x is the eigenvalue.
%   EIGSHOW, by default, uses
%      A = [1/4 3/4; 1 1/2];
%   Other interesting examples are:
%      A = eye(2,2);
%      A = [0 1; 1 0];
%      A = [0 1; -1 0];
%      A = [1/2 1; -1/4 -1/2]
%      A = [3/2 1; -1/4 1/2]
%      A = randn(2,2);
%
%   See also SVDSHOW.

%   Copyright (c) 1993-97 by The MathWorks, Inc.
%   $Revision: $  $Date: $

if nargin == 0;
   A = [1/4 3/4; 1 1/2];
   eigshowinit(A);
elseif all(size(A) == 2)
   eigshowinit(A);
elseif isequal(A,'action')
   eigshowaction
else
   error('Matrix must be 2-by-2.')
end

%------------------

function eigshowinit(A)
clf reset
s = 1.1*max(1,norm(A));
axis([-s s -s s])
axis square
title(sprintf('A = [ %6.3f  %6.3f; %6.3f  %6.3f ]',A'))
xlabel('Make A*x parallel to x')
xcolor = [0 .6 0];
Axcolor = [0 0 .8];
h.x = eigshowvec([1 0]','x',xcolor);
h.Ax = eigshowvec(A(:,1),'Ax',Axcolor);
set(gcf,'userdata',A);
set(gca,'userdata',h);
set(gcf,'windowbuttondownfcn', ...
  'eigshow(''action''); set(gcf,''windowbuttonmotionfcn'',''eigshow(''''action'''')'')')
set(gcf,'windowbuttonupfcn', ...
  'set(gcf,''windowbuttonmotionfcn'','''')')

%------------------

function eigshowaction
A = get(gcf,'userdata');
h = get(gca,'userdata');
pt = get(gca,'currentpoint');
x = pt(1,1:2)';
x = x/norm(x);
eigshowact(x,h.x);
eigshowact(A*x,h.Ax);

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
