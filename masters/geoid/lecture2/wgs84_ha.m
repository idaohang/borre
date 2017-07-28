function N = wgs84_ha(phi,lambda,H,n_max)
%WGS84_HA Computes the geoidal undulation N from a spherical harmonic
%         expansion of the Earth's gravity field for given geographical
%         coordinates (phi,lambda) in degrees and orthometric height
%         H in meters. The series expansion is summed up to and 
%         including n_max.
%         We use the EGM96 spherical harmonic coefficients as given
%         in the file egm180.nor. (We have added a line for n=2, m=1.)
%         This file can be found at 
%         64.214.2.59/GandG/wgs-84/user-cle.html

% References: Heiskanen, Weikko A. & Helmut Moritz (1967):
%             Physical Geodesy. W.H.Freeman and Company
%
%             Department of Defence World Geodetic System 1984.
%             Third edition, 1997. National Imagery and Mapping 
%             Agency. Technical Report

%Written by Kai Borre
%March 25, 2000

if n_max > 180, disp('Upper summation bound too large!'), end
if nargin == 0
   phi = 57;
   lambda = 11;
   H = 0;
   n_max = 18;
end

GM = 3986004.418e8;  % m^3/s^2
a = 6378136.46;      % m
%we input lambda and phi in degrees:
[X,Y,Z] = frgeod(6378137,298.257223563,phi,lambda,H);
[lambda,phi_,r] = cart2sph(X,Y,Z); % lambda and phi are now in radians
sin_phi_ = sin(phi_); % phi_ is geocentric latitude
sn = 0;
fid = fopen('egm180.nor','r');
for n = 2:n_max
   sm = 0;
   for m = 0:n
      line = fgets(fid);
      n_out = line(1:5);
      m_out = line(6:10);
      % fprintf('\n n_out %5.0f m_out %5.0f n %5.0f',... 
      %                      str2num(n_out),str2num(m_out),n)
      Cnm = str2num(line(11:25));
      Snm = str2num(line(26:40));
      if m == 0
         switch n
         case 2,  Cnm = Cnm+0.484166774985e-3;
         case 4,  Cnm = Cnm-0.790303733511e-6;
         case 6,  Cnm = Cnm+0.168724961151e-8;
         case 8,  Cnm = Cnm-0.346052468394e-11;
         case 10, Cnm = Cnm+0.265002225747e-14; 
         end  
      end     
      sm = sm+(Cnm*cos(m*lambda)+Snm*sin(m*lambda))* ...
                                   legen_nm(n,m,sin_phi_);
   end
   sn = sn+(a/r)^n*sm;   
end
phi = phi*pi/180;
zeta = GM*sn/(r*gamma_h(phi,H));  % height anomaly
if H ~= 0
   h = H+zeta;
   zeta = GM*sn/(r*gamma_h(phi,h));
end
N = -0.53+zeta;
fclose(fid);
fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%% end wgs84_ha %%%%%%%%%%%%%%%