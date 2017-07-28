function u = fmintx(F,a,b,tol,varargin)
%FMINTX  Textbook version of FMINBND
%   x = FMINTX(F,a,b) finds a local minimizer x of the function F
%   in the interval a <= x <= b. F accepts scalar input x and returns 
%   a scalar function value, F(x).
%
%   x = FMINTX(F,a,b,tol) uses stopping tolerance tol instead of 1.e-6.
%
%   x = FMINTX(F,a,b,tol,p1,p2,...) provides for additional
%   arguments, which are passed to the objective function, F(x,p1,p2,...).
%   (Use [] as a place holder for tol to get the default tolerance.)
%
%   Examples
%     F can be specified using @:
%        x = fmintx(@cos,3,4)
%      computes pi to six decimal places.
%
%     F can also be an inline object:
%        f = inline('sin(x)+3');
%        x = fmintx(f,2,5);
%
%   See also FMINBND, FMINSEARCH, FZERO, @, INLINE.

%   Reference: "Computer Methods for Mathematical Computations",
%   Forsythe, Malcolm, and Moler, Prentice-Hall, 1976.

% Initialization

if nargin < 4 | isempty(tol)
   tol = 1.e-6; 
end
if ischar(F) & exist(F)~=2
   F = inline(F);
elseif isa(F,'sym')
   F = inline(char(F));
end 

phi = (1 + sqrt(5))/2;
rho = 2 - phi;
u = a + rho*(b-a);
v = u; w = u; x = u;
fu = feval(F,u,varargin{:}); 
fv = fu; fw = fu; fx = fu;
xm = 0.5*(a+b);
d = 0.0;
e = 0.0;

% Main loop

while abs(x-xm) > tol
   % Is a parabolic fit possible?
   para = abs(e) > tol;
   if para
      % Try parabolic fit.
      r = (x-w)*(fx-fv);
      q = (x-v)*(fx-fw);
      p = (x-v)*q-(x-w)*r;
      q = 2.0*(q-r);
      if q > 0.0,  p = -p; end
      q = abs(q);
      r = e;
      e = d;
      % Is the parabola acceptable?
      para = ( (abs(p)<abs(0.5*q*r)) & (p>q*(a-x)) & (p<q*(b-x)) );
      if para
         % Parabolic interpolation step
         d = p/q;
      end
   end
   if ~para
      % Golden-section step
      if x >= xm
         e = a-x;
      else
         e = b-x;
      end
      d = rho*e;
   end
   
   u = x + d;
   fu = feval(F,u,varargin{:});  
   
   % Update a, b, x, v, w, xm
   if fu <= fx
      if u >= x, a = x; else, b = x; end
      v = w; fv = fw;
      w = x; fw = fx;
      x = u; fx = fu;
   else
      if u < x, a = u; else, b = u; end
      if ( (fu <= fw) | (w == x) )
         v = w; fv = fw;
         w = u; fw = fu;
      elseif ( (fu <= fv) | (v == x) | (v == w) )
         v = u; fv = fu;
      end
   end
   xm = 0.5*(a+b);
end
