%FFT_GEOI  Computation of local geoid in Danmark from free-air
%          gravity anomalies at more than 40 000 points

% Written by Kai Borre
% April 19, 2000

%fid = fopen('dk.dat','r')
%D = fscanf(fid,'%g %g %g %g %g',[5 inf]); 
load D
D = D';
[m,n] = size(D);
%Phi = D(:,2);
%Lambda = D(:,3);
%Delta_g = D(:,5);
%%figure(1);
%%plot(D(:,3),D(:,2),'.')
%%title(['Distribution of ' int2str(m) ' Danish gravity data'])

ml = 16; % x-direction
nl = 16; % y-direction
ylin = linspace(55,56,nl); % we interpolate between phi = 56 and
                           % 57, lambda = 9 and 10
xlin = linspace(9,13,ml);
Deltax = 111000*cos(56.5*pi/180)/ml; % meter
Deltay = 111000/nl;                 % meter
[X,Y] = meshgrid(xlin,ylin);
Z = griddata(D(:,3),D(:,2),D(:,5),X,Y);     % mgal

%figure(2);
%mesh(X,Y,Z);
%hold on
%plot3(D(:,3),D(:,2),D(:,5),'.')
%axis([9 10 56 57 -50 50])
%hold off

tic
N = zeros(ml,nl);
FG = fft2(Z);
for k = 1:ml
   k
   for l = 1:nl
      L = distance(k,l,ml,nl,Deltax,Deltay);
      FL = fft2(L);
      FP = FG*FL;
      geo = ifft2(FP);
      deltag = sqrt(Deltax*Deltay)*Z(k,l)/(sqrt(pi)*9.8);
      N(k,l) = geo(k,l)+ deltag;  
   end
end
N = N/(2*pi*9.8);
figure(2);
[C,h] = contour(N);
clabel(C,h);
toc
%%%%%%%%%%%%%%%%%%%%% end fft_geoi.m  %%%%%%%%%%%%%%%%%%