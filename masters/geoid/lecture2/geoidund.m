function N = geoidund(phi,lambda,H,n_max)
%GEOIDUND Computes the geoidal undulation N from a spherical harmonic
%         expansion of the Earth's gravity field for given geographical
%         coordinates (phi,lambda) in degrees and orthometric height
%         H in meters. The series expansion is summed up to and 
%         including n_max.
%         We use the EGM96 spherical harmonic coefficients as given
%         in the file egm180.nor. (We have added a line for n=2, m=1.)
%         This file can be found at 
%         www.nima.mil/GandG/wgsegm/egm96.html

% References and recommended reading:
%             Heiskanen, Weikko A. & Helmut Moritz (1967):
%             Physical Geodesy. W.H.Freeman and Company
%
%             Department of Defence World Geodetic System 1984 (1997):
%             Third edition. National Imagery and Mapping Agency. 
%             Technical Report
%
%             Holmes, S.A. & W.E. Featherstone (2002): A unified
%             approach to the Clenshaw summation and the recursive
%             computation of very high degree and order normalised 
%             associated Legendre functions. Journal of Geodesy 76: 
%             279--299
%        
%             Bj\"o{}rck, \AA{}ke (1996): Numerical Methods for Least 
%             Squares Problems. SIAM

% Based on NIMA's public code of 10 October 1996:
% earth-info.nima.mil/GandG/wgsegm/clenqt.for

% Written by Kai Borre, April 15, 2004

if nargin == 0
    phi = 57;
    lambda = 10;
    H = 0;
    n_max = 120;
end
if n_max > 180, disp('Upper summation bound too large!'), return, end

GM = 3986004.418e8;     % m^3/s^2
a_e = 6378137.0;        % m
finv = 298.257223563;   % []

% We input phi and lambda in degrees
[X,Y,Z] = frgeod(a_e,finv,phi,lambda,H);
[lambda,phi_geoc,r] = cart2sph(X,Y,Z); % lambda and phi_geoc are in radians
t = sin(phi_geoc); % phi_geoc is geocentric latitude
y = cos(phi_geoc);

%We prepare the double summation by choosing a column oriented 
%sequence of the terms. The actual summation of terms will be 
%performed as a Clenshaw summation, see Bj{\"o}rck (1996) eq. (8.2.7)

if exist('geoc.dat') ~= 2
    Cnm = zeros((n_max+1)*(n_max+2)/2,1);
    Snm = zeros((n_max+1)*(n_max+2)/2,1);
    % Reading the Normalized Geopotential Coefficients Cnm and Snm
    disp('Reading Coefficient File')
    fid = fopen('egm180.nor','r');
    for n = 2:n_max
        for m = 0:n   
            n = fscanf(fid,'%d',1);
            m = fscanf(fid,'%d',1);
            nm = n*(n+1)/2+m+1;
            c = fscanf(fid,'%15cs',1);
            Cnm(nm) = str2double(c);
            s = fscanf(fid,'%15cs',1);
            Snm(nm) =str2double(s);
            %fprintf('\n ord %d deg %d:   Cnm %-9.8E   Snm %-9.8E', ...
            %             n,m,Cnm(nm),Snm(nm))
        end
    end
    fclose(fid);
    %Subtracting the normal gravity field by
    % modifying even zonal coefficients Cnm.          
    % WGS84 values are used
    esq = 0.00669437999013; % []
    C2 = 108262.9989050e-8; % []
    for i = 2:5
        C2n(i,1) = (-1)^(i+1)*(3*(esq^i)/((2*i+1)*(2*i+3)))*(1-i+5*i*C2/esq);
    end
    Cnm(4) = Cnm(4)+C2/sqrt(5);
    Cnm(11) = Cnm(11)+C2n(2)/3;
    Cnm(22) = Cnm(22)+C2n(3)/sqrt(13);
    Cnm(37) = Cnm(37)+C2n(4)/sqrt(17); 
    Cnm(56) = Cnm(56)+C2n(5)/sqrt(21);
    save geoc.dat Cnm -ascii -double 
    save geos.dat Snm -ascii -double
else
    Cnm = load('geoc.dat');
    Snm = load('geos.dat');
end
disp('Completed Reading Coefficients')

% Computing sin(j*lambda) and cos(j*lambda) by recursion from 
% powers of sin(lambda) and cos(lambda)
sinm = zeros(n_max+1,1); 
cosm = zeros(n_max+1,1);
% sinm(1) = 0;
cosm(1) = 1;
sinm(2) = sin(lambda);
cosm(2) = cos(lambda);
for j = 3:n_max+1
    sinm(j) = 2*cosm(2)*sinm(j-1)-sinm(j-2); % sin(j*lambda)
    cosm(j) = 2*cosm(2)*cosm(j-1)-cosm(j-2); % cos(j*lambda)
end

% build all Clenshaw coefficient arrays
disp('Building Clenshaw Coefficient Arrays')
as = zeros(n_max+3,1); 
for i = 1:n_max+1
    as(i) = -sqrt((2*i+1)/(2*i));
end    

% building location array
loc = zeros(n_max+3,1);
for i = 1:n_max+3
    loc(i)=i*(i-1)/2+1;
end
% Normalisation coefficients a_nm and b_nm
a = zeros((n_max+3)*(n_max+4)/2,1); 
b = a;
for m = 1:n_max+1 
    for n = m+1:n_max 
        s = loc(n+1)+m; 
        a(s) = sqrt(((2*n+1)*(2*n-1))/((n-m)*(n+m)));
        b(s) = sqrt(((2*n+1)*(n+m-1)*(n-m-1))/((2*n-3)*(n+m)*(n-m)));
    end
end

q = a_e/r;
s1 = zeros(n_max+4,1); 
s2 = s1;
sht = s1;
for m = n_max:-1:2 
    for n = n_max:-1:m
        loc1 = loc(n+1)+m; 
        loc2 = loc(n+2)+m; 
        loc3 = loc(n+3)+m; 
        s1(n) = a(loc2)*t*q*s1(n+1) - b(loc3)*q^2*s1(n+2) + Cnm(loc1); 
        s2(n) = a(loc2)*t*q*s2(n+1) - b(loc3)*q^2*s2(n+2) + Snm(loc1);        
    end
    sht(m) = -as(m+1)*y*q*sht(m+1) + s1(m)*cosm(m) + s2(m)*sinm(m);
end
% The final result from the Clenshaw summation
sn = (s1(1)+s2(1))*q + sht(2)*sqrt(3)*y*q^2; 
phi = phi*pi/180;
zeta = GM*sn/(r*gamma_h(phi,H));  % height anomaly
if H ~= 0
    h = H+zeta;
    zeta = GM*sn/(r*gamma_h(phi,h)); 
end
N = zeta-0.53;
fprintf('\n Geoidal Undulation N = %6.3f meters\n\n',N);
%%%%%%%%%%%%%%%%%%%%%%% end geoidund %%%%%%%%%%%%%%%