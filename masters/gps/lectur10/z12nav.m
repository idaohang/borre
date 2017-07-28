function rawnav = z12nav(fidobs)
%Z12NAV   Reads binary observations in an already opened and positioned 
%         file; the observations are from a Z-12 receiver.

%Kai Borre 17-03-03
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2003/03/17 $

% Initial computations of constants
v_light = 299792458; 	     % vacuum speed of light m/s
% The size of rawnav structure = 67 bytes        
RA = fread(fidobs,4,'char');
RAs = setstr(RA');
if feof(fidobs) == 1, return, end
RB = fread(fidobs,4,'double'); % 8char
RC = fread(fidobs,3,'float');  % 4char
RD = fread(fidobs,2,'double'); % 8char
RE = fread(fidobs,1,'uint16'); % 2char
RF = fread(fidobs,1,'char');
%Each epoch has a rawdata struct per SV
%The size of L1  rawdata struct = (35 * nav.num_sats) bytes
%The size of L2C rawdata struct = (66 * nav.num_sats) bytes
%The size of L2P rawdata struct = (97 * nav.num_sats) bytes
% The proper observations:
% The counter j runs through the individual satellites in each epoch
% The counter i runs through the various code observations
% the number of which is rawheader.num_obs_types. 
Asv = zeros(32,1)*nan;
Ael = zeros(32,1)*nan;
Aaz = zeros(32,1)*nan;
Ach = zeros(32,1)*nan;
Ai = zeros(32,3)*nan;
Bi = zeros(32,3)*nan;
Ci = zeros(32,3)*nan;
Di = zeros(32,3)*nan;
E1i = zeros(32,3)*nan;
E2i = zeros(32,3)*nan;
E3i = zeros(32,3)*nan;
Fi = zeros(32,3)*nan;
Gi = zeros(32,3)*nan;
Hi = zeros(32,3)*nan;    
for j = 1:RF
    % i = 1:   C/A-code on L1
    % i = 2: 	 P-code on L1
    % i = 3:     P-code on L2
    s1 = fread(fidobs,1,'uchar');
    sv = s1;
    Asv(sv,1) = s1;
    Ael(sv,1) = fread(fidobs,1,'uchar');
    Aaz(sv,1) = fread(fidobs,1,'uchar');
    Ach(sv,1) = fread(fidobs,1,'uchar');
    for i = 1:3
        Ai(sv,i) = fread(fidobs,1,'double');
        Bi(sv,i) = fread(fidobs,1,'float');
        Ci(sv,i) = fread(fidobs,1,'ushort');
        Di(sv,i) = fread(fidobs,1,'char');
        E1i(sv,i) = fread(fidobs,1,'uchar');  
        E2i(sv,i) = fread(fidobs,1,'uchar');
        E3i(sv,i) = fread(fidobs,1,'uchar');
        Fi(sv,i) = fread(fidobs,1,'char');
        Gi(sv,i) = fread(fidobs,1,'long');
        Hi(sv,i) = fread(fidobs,1,'double');    
    end % i
end % j
rawnav = struct('sitename',RAs(1,1:4), ...
    'rcv_time',RB(1,1), ...
    'navx',RB(2,1), ...
    'navy',RB(3,1), ...
    'navz',RB(4,1), ...
    'navxdot',RC(1), ...
    'navydot',RC(2), ...
    'navzdot',RC(3), ...
    'navt',RD(1), ...
    'navtdot',RD(2), ...
    'pdop',RE(1)/100, ...
    'num_sats',RF(1), ...
    'nest', ...
    struct('svprn',Asv(:,1),...
    'elevation',Ael(:,1),...
    'azimuth',2*Aaz(:,1),...
    'chnind',Ach(:,1),...       
    'rawrange',Ai(:,1)*v_light,... % SV raw range: raw transmit time is this  
    ...                       % value subtracted from receive time       
    ...                       % (receive time rounded to nearest         
    ...                       % millisecond)       
    'smth_corr',Bi(:,1),...        % Smoothing correction for ranges (meters)
    'smth_count',Ci(:,1),...       % Number of data points in smoothing   
    'polarity_known',Di(:,1),...   % RESERVED   
    'warning',E1i(:,1),...         % Warning flag (BIT flags):                 
    ...                       %  Bit 1 ==> RESERVED                       
    ...                       %  Bit 2 ==> RESERVED                       
    ...                       %  Bit 3 ==> Carrier phase questionable.    
    ...                       %  Bit 4 ==> Code phase questionable.       
    ...                       %  Bit 5 ==> RESERVED                       
    ...                       %  Bit 6 ==> PW tracking method used.       
    ...                       %  Bit 7 ==> Possilbe loss of lock.         
    ...                       %  Bit 8 ==> Loss of lock occured.          
    'goodbad',E2i(:,1),...         % Another health indicator:                
    ...                       %  22 ==> Code and carrier measured.        
    ...                       %  23 ==> Same as 22 but additionally, nav  
    ...                       %         message obtained but measurement 
    ...                       %         was not used in position         
    ...                       %         computation.                     
    ...                       %  24 ==> Same as 23 but codephase was used 
    ...                       %         in position computation.         
    'ireg',E3i(:,1),...            % Signal to noise
    'qa_phase',Fi(:,1),....        % QA phase check (0.001 cycles)
    'doppler',Gi(:,1),...          % SV raw doppler
    'carphase',Hi(:,1)));           % Full carrier phase (cycles)    
for q = 2:3
    rawnav.nest(q).rawrange = Ai(:,q)*v_light;      
    rawnav.nest(q).smth_corr = Bi(:,q);
    rawnav.nest(q).smth_count = Ci(:,q);   
    rawnav.nest(q).polarity_known = Di(:,q);   
    rawnav.nest(q).warning = E1i(:,q);                 
    rawnav.nest(q).goodbad = E2i(:,q);       
    rawnav.nest(q).ireg = E3i(:,q);
    rawnav.nest(q).qa_phase = Fi(:,q);
    rawnav.nest(q).doppler = Gi(:,q);
    rawnav.nest(q).carphase = Hi(:,q);
end % q
%%%%%%%% end z12nav.m %%%%%%%%%%%%%%%%%%%%%%%%

