% CLENCT  A transcribed Fortran code (NIMA/GIMG, 10 OCT 1996)        
%         for computing quantities related to the 
%         WGS84 EM96 geoid. For numerical reasons the summation is 
%         done according to Clenshaw's recursion formula,
%         cf. Aake Bjoerck (1996): Numerical Method for Least 
%         Squares Problems, SIAM.
%         NMAX denotes maximum degree and order attained.
%         ND is number of C or S geopotential coefficients.
%         NDD is size of A and B Clenshaw arrays.
%         We use the EGM96 spherical harmonic coefficients as given
%         in the file egm180.nor. (We have added a line for n=2, m=1.)
%         This file can be found at 
%         http://earth-info.nga.mil/GandG/wgs84/gravitymod/wgs84_180/user-cle.html

% References:  Department of Defence World Geodetic System 1984.
%             Third edition, 1997. National Imagery and Mapping 
%             Agency. Technical Report

% Fortran code (CLENCT.FOR) transcribed by Kai Borre
% January 2, 2001

% Revised by Kostas Dragunas
% September 25, 2006


% Example of call
%  out = clenct([45.5 35.5 -30], [15 10 310.1 ], [0 2 5]);
%  [out in] = clenct([45.5 35.5 -30], [15 10 310.1 ], [0 2 5]);

function [out in] = clenct(PHID, RLAMD, HT)

% input:
%   PHI      - Latitude in degrees in the range (-90 90).
%              The poles are not included
%              Multiple calls can be arranged as a column or a 
%              row vector
%   LAMBDA   - Longitude in degrees in the range [0 360]
%              Multiple calls can be arranged as a column or 
%              a row vector
%   HT       - Geodetic Height in meters. Multiple calls can 
%              be arranged as a column or a row vector
%              
%              The code checks the dimensions of PHI, LAMDA
%              and HT and will generate output according to this
%              information. User inputs the variables PHID, RLAMD, 
%              and HT of the same dimensions (e.g. 1xN, 1xN, 1xN). 
%              There is no difference if input is arranged as a 
%              column or a row vector (e.g. will 1xN, Nx1, 1xN 
%              work). If the size of input parameters are
%              different, the output will be as big as:
%              out_size = min([length(PHID) length(RLAMD) length(HT)])              
% 
% output:
%   out      - matrix with the following parameters:
%               SUMHT - point geoid height, N in meters, 
%               SUMG  - point free-air gravity anomaly, DELTA G in mgal,
%               DR    - point radial component of the gravity disturbance vector 
%                       DIST R in mgal,
%               DNS   - point north-south components of the gravity
%                       disturbance vector, DIST NS in mgal,
%               DEW   - point east-west components of the gravity disturbance vector 
%                       DIST EW in mgal,
%               SUMX  - point north-south components of the deflection of
%                       the vertical, XI in arcsec, 
%               SUME  - point east-west components of the deflection of
%                       the vertical, ETA in arcsec,
%               DEFL  - point total deflection of the vertical, THETA in
%                       arcsec.
%              Each row in matrix 'out' contains results for
%              different point. 

%   in       - matrix with the input parameters. This parameter is
%              optional, however, if the input parameters are of different 
%              length, the user will not know exactly which input was used 
%              to compute which points. Each row of this matrix represents
%              the point parameters PHID, RLAMD, and HT. Each row in "out" 
%              matrix matches the same row in the "in" matrix.
%             
%
% This script uses the data file EGM180.NOR (or optimized EGM180_opt.nor). 
% It contains the WGS 84 n = m = 180, normalized, geopotential coefficients.       
% The file EGM180.NOR is an ASCII file rowwise containing the degree (n), order (m), 
% and normalized geopotential coefficients (Cnm, Snm) in the following format: 
%                n m Cnm Snm 
% These coefficients were developed using surface gravity, satellite altimetry, and 
% elevation data. 
% More information about the file:
% http://earth-info.nga.mil/GandG/wgs84/gravitymod/wgs84_180/user-cle.html

%%global data_in

% default values for input and output
out = [];
in  = [];

% for time computation
time_x = clock;

% allocation arrays
NMAX = 180;
NP1  = NMAX+1;
NP2  = NMAX+2;
ND   = NP1*NP2/2;
NP3  = NMAX+3;
NP4  = NMAX+4;
NDD  = NP3*NP4/2;

% OTHER ARRAYS STORE CLENSHAW COEFFICIENTS NEEDED FOR
% THE SELECTED GRAVIMETRIC QUANTITIES THAT ARE COMPUTED
T11  = zeros(1,NP2+1);
TG1  = zeros(1,NP2+1);
TP1  = zeros(1,NP2+1);
T12  = zeros(1,NP2+1);
TG2  = zeros(1,NP2+1);
TP2  = zeros(1,NP2+1);
S11  = zeros(1,NP2+1);
S12  = zeros(1,NP2+1);
SP1  = zeros(1,NP2+1);
SP2  = zeros(1,NP2+1);
SV15  = zeros(1,NP2+1);
SV1  = zeros(1,NP2+1);
SV2  = zeros(1,NP2+1);
RNN  = zeros(1,NMAX+1);
A    = zeros(1,NDD); 
B    = zeros(1,NDD);
CNM  = zeros(1,ND); % store geopotential coefficients
SNM  = zeros(1,ND); % store geopotential coefficients
IV   = zeros(1,NP3); % locating array
AS   = zeros(1,NP3+1);
CR   = zeros(1,NP1+1);
SR   = zeros(1,NP1+1);
S2   = zeros(1,NP3+1);
SG1  = zeros(1,NP3+1);
SG2  = zeros(1,NP3+1);
SG   = zeros(1,NP3+1);
SHT  = zeros(1,NP3+1); 
C2N  = zeros(1,5); 

% WGS84 MODEL VALUES ARE USED
C2    = 108262.9989050e-8;
RKM   = 3.98600150e14;
AE    = 6378137.0;
ESQ   = 0.006694379990130;
%F     = 298.257;
GRAVA = 9.7803267714;
STAR  = 0.001931851386;

% if WGS72 MODEL VALUES ARE USED
%   C2    = 108263.0e-8;
%   RKM   = 3.9860050e14;
%   AE    = 6378135.0;
%   ESQ   = 0.0066943177780;
%   F     = 298.260;
%   GRAVA = 9.78033270;
%   STAR  = 0.0052789940;

TEMPP  = 400.0;
TEMPL  = 400.0;
HGHT   = 0.0000001;
GAMMAA = RKM/AE^2;
GA     = RKM/AE;  
T1     = 1.0;
T2     = 2.0;
T3     = 3.0;
PI2    = pi/2.0;
RPD    = pi/180.0;  
RPS    = RPD/3600.0;
RHO    = 1.0/RPS;

% READING POTENTIAL COEFFICIENTS
% Reading the modified data file which is a space separated file.
% The reason for reformatting the original file is twofold:
%  * the file becomess shorter
%  * the file is faster to read (~6x)
% 
% If this type of file can not be used for some reason, then comment 
% the following 22 lines and uncomment the code below which will do 
% the same thing as the original Fortran code

%%if exist('data_in','var') == 0
    
%fprintf('\n- Reading Coefficient File (egm180_opt.nor)... ')
fid = fopen('egm180_opt.nor','r');
if fid == -1
    fprintf('\n- File "egm180_opt.nor" have not been found!');
    return;
end
N = 0;
tic
while N <= NMAX
    line = fgets(fid);
    if line ~= -1
        data_in = sscanf(line,'%d %d %f %f')';
        N  = data_in(1);
        M  = data_in(2);
        LL = N * (N + 1) / 2 + M + 1;
        CNM(LL) = data_in(3);
        SNM(LL) = data_in(4);
    else
        break;
    end
end
fclose(fid);
%fprintf(['\n- Finished Reading Coefficients '...
%         'in %.2f seconds'],toc);
%end

% File reading procedure from the old data type file.
% This is a slow procedure so it is better to use the 
% one shown above
%
% fprintf('\n- Reading Coefficient File (egm180.nor)... ')
% fid = fopen('egm180.nor','r');
% if fid == -1
%     fprintf('\n- File "egm180.nor" has not been found!');
%     return;
% end
% N = 0;
% tic
% while N <= NMAX
%     line = fgets(fid);
%     if line ~= -1
%         N  = str2double(line(1:5));
%         M  = str2double(line(6:10));
%         LL = N * (N + 1) / 2 + M + 1;
%  
%%%% CNM(LL) = str2double(line(11:25));
%         SNM(LL) = str2double(line(26:40));
%     else
%         break;
%     end
% end
% fclose(fid);

SR(1) = 0;
CR(1) = 1;

% BUILD LOCATING ARRAY
for I = 1:NP3
   IV(I) = I * (I - 1) / 2;
end 
for I = 0:NMAX 
   RNN(I+1) = I - 1;
end 

% MODIFY CNM EVEN ZONAL COEFFICIENTS  
for I = 2:5  
   C2N(I) = (-1.0)^(I+1) * (3.0 * (ESQ^I) / ((2.0 * I + 1) * (2.0 * I + 3))) * (1 - I + 5.0 * I * C2 / ESQ);
end 

% WGS84 VALUES ARE USED
CB2     = -1 * C2 / sqrt(5.0); 
CB4     = -1 * C2N(2) / 3.0;
CB6     = -1 * C2N(3) / sqrt(13.0);
CB8     = -1 * C2N(4) / sqrt(17.0);
CB10    = -1 * C2N(5) / sqrt(21.0);

% if WGS72 VALUES ARE USED. 
%      CB2 = -4.841732E-04;  
%      CB4 = 7.8305E-07; 
%      CB6 = 0.0;
%      CB8 = 0.0;
%      CB10 = 0.0;
  
CNM(4)  = CNM(4) - CB2;
CNM(11) = CNM(11) - CB4;  
CNM(22) = CNM(22) - CB6;  
if NMAX > 6, CNM(37) = CNM(37) - CB8; end
if NMAX > 9, CNM(56) = CNM(56) - CB10; end

% BUILD ALL CLENSHAW COEFFICIENT ARRAYS
for I = 1:NMAX
   RI = I;
   AS(I+1) = -1.0 * sqrt((T2 * RI + T1) / (T2 *RI));  
end 
for M = 0:NMAX  
   MP1 = M+1;  
   RM = M;
   for N = MP1:NMAX
      RN = N;
      LL = IV(N+1) + M + 1; 
      A(LL) = -1.0 * sqrt(((T2 * RN + T1) * (T2 * RN - T1)) / ((RN - RM) * (RN + RM)));
      B(LL) = sqrt(((T2 * RN + T1) * (RN + RM - T1) * (RN - RM - T1)) / ((T2 * RN - T3) * (RN + RM) * (RN - RM)));
   end 
end 

size_input = min([length(PHID) length(RLAMD) length(HT)]);

% initialize output
out = zeros(size_input,8);

% initialize input
in = zeros(size_input,3);

for ii = 1:size_input
    in(ii,:) = [PHID(ii) RLAMD(ii) HT(ii)];
    if abs(PHID(ii)) < 90 &&  RLAMD(ii)>=0 && RLAMD(ii)<=360

        PHI   = PHID(ii)*RPD;

        % PERFORM CLENSHAW SUMMATION FOR EACH PT
        % COMPUTE GEOCENTRIC DISTANCE, COLATITUDE AND NORMAL GRAVITY
        if (TEMPP ~= PHID(ii)) && (HGHT ~= HT(ii))
            RN   = AE / sqrt(1 - ESQ * (sin(PHI))^2);
            T22  = (RN + HT(ii)) * cos(PHI);
            X2Y2 = T22^2;
            Z1   = (RN * (1 - ESQ) + HT(ii)) * sin(PHI);
            R    = sqrt(X2Y2 + Z1^2);
            PSI  = atan2(Z1,sqrt(X2Y2));

            %WGS84 VALUES ARE USED
            GRAVN = GRAVA * (1 + STAR * (sin(PHI))^2) / sqrt(1 - ESQ * (sin(PHI))^2);
            % if WGS72 VALUES ARE USED.
            %    GRAVN = GRAVA * (1 + STAR * DSIN(PHI)^2 + 0.000023461 * DSIN(PHI)^4);

            GRAVN = GRAVN - HT(ii) * 0.3086e-5;
            DC    = cos(PSI);
            TH    = PI2-PSI;
            Y     = sin(TH);
            T     = cos(TH);
            XX    = T/Y;
            F1    = AE/R;
            F2    = (AE/R)^2;
            P11   = sqrt(T3) * Y * F1^3;
            P11HT = sqrt(T3) * Y * F2;
            D11   = -1.0 * sqrt(T3) * XX * F1^3;
        end

        if TEMPL ~= RLAMD(ii)
            RLAM = RLAMD(ii) * RPD;
            SR(2) = sin(RLAM);
            CR(2) = cos(RLAM);
            for J = 2:NMAX
                SR(J+1) = T2 * CR(2) * SR(J) - SR(J-1);
                CR(J+1) = T2 * CR(2) * CR(J) - CR(J-1);
            end
        end

        for M = NMAX:-1:0
            if TEMPP ~= PHID(ii) && HGHT ~= HT(ii)
                for N = NMAX:-1:M
                    LL  = IV(N+1) + M + 1;
                    LL2 = IV(N+2) + M + 1;
                    LL3 = IV(N+3) + M + 1;
                    S11(N+1) = -1.0 * A(LL2)*T*F1*S11(N+2)-B(LL3)*F2*S11(N+3)+CNM(LL);
                    S12(N+1) = -1.0 * A(LL2)*T*F1*S12(N+2)-B(LL3)*F2*S12(N+3)+SNM(LL);
                    SG1(N+1) = -1.0 * A(LL2)*T*F1*SG1(N+2)-B(LL3)*F2*SG1(N+3)+CNM(LL)*RNN(N+1);
                    SG2(N+1) = -1.0 * A(LL2)*T*F1*SG2(N+2)-B(LL3)*F2*SG2(N+3)+SNM(LL)*RNN(N+1);
                    SP1(N+1) = -1.0 * A(LL2)*T*F1*SP1(N+2)-B(LL3)*F2*SP1(N+3)-A(LL2)*F1* S11(N+2);
                    SP2(N+1) = -1.0 * A(LL2)*T*F1*SP2(N+2)-B(LL3)*F2*SP2(N+3)-A(LL2)*F1* S12(N+2);
                end
                T11(M+1) = S11(M+1);
                T12(M+1) = S12(M+1);
                TG1(M+1) = SG1(M+1);
                TG2(M+1) = SG2(M+1);
                TP1(M+1) = SP1(M+1);
                TP2(M+1) = SP2(M+1);
            end

            SV15(M+1) = -1.0 * AS(M+2)*Y*F1*SV15(M+2)+T11(M+1)*CR(M+1)+T12(M+1)*SR(M+1);
            SV1(M+1)  = -1.0 * AS(M+2)*Y*F1* SV1(M+2)+TP1(M+1)*CR(M+1)+TP2(M+1)*SR(M+1);
            SV2(M+1)  = -1.0 * AS(M+2)*Y*F1* SV2(M+2)+ AS(M+2)*XX*F1*SV15(M+2);
            S2(M+1)   = -1.0 * AS(M+2)*Y*F1*  S2(M+2)-T11(M+1)*M*SR(M+1)+T12(M+1)*M*CR(M+1);
            SG(M+1)   = -1.0 * AS(M+2)*Y*F1*  SG(M+2)+TG1(M+1)*CR(M+1)+TG2(M+1)*SR(M+1);
            SHT(M+1)  = -1.0 * AS(M+2)*Y*F1* SHT(M+2)+T11(M+1)*CR(M+1)+T12(M+1)*SR(M+1);
        end

        SUMX  = (SV1(2)+SV2(2))*P11 + SV15(2)*D11 +(TP1(1)+TP2(1))*F2;
        SUMX  = (-1.0*Y*SUMX*GAMMAA/GRAVN)*RHO;
        SUME  = S2(2)*P11;
        SUME  = -(GAMMAA*SUME*RHO)/(GRAVN*DC);
        SUMG  = (TG1(1)+TG2(1))*F2+SG(2)*P11;
        SUMG  = SUMG*GAMMAA*1.0E5;
        SUMHT = (T11(1)+T12(1))*F1+SHT(2)*P11HT;
        SUMHT = SUMHT*GA/GRAVN;
        DR    = -1.0 * (SUMG+2.0*SUMHT*GRAVN/R*1.0E5);
        DNS   = -1.0 * GRAVN*SUMX/RHO*1.0E5;
        DEW   = -1.0 * GRAVN*SUME/RHO*1.0E5;
        DEFL  = sqrt(SUMX^2+SUME^2);
        
        out(ii,:) = [SUMHT, SUMG, DR, DNS, DEW, SUMX, SUME, DEFL];

     else
        out(ii,:) = [0, 0, 0, 0, 0, 0, 0, 0];
        if abs(PHID(ii)) >= 90 
            fprintf(['\n- You are close to a Pole,' ...
                     'thus, no output for this point' ....
                      'PHID(%d)= %f should be (-90 90)'], ii,PHID(ii));
        elseif RLAMD(ii)<0 || RLAMD(ii)>360
            fprintf(['\n- Longitude exceeds the limits,'...
                     'thus, no output for this point.'...
                     'RLAMD(%d)= %f should be [0 360]'], ii,RLAMD(ii));
        end
    end
end
    
%fprintf(['\n- Program finished computations in ' ...
%         '%.2f seconds (total time)\n'], etime(clock,time_x));

%%%%%%%%%%%%%%%%%%%%%%% end clenct.m %%%%%%%%%%%