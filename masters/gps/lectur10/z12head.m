function [rawheader,fidobs] = z12head(observationfile)
%Z12HEAD  Reads the header of a binary data file created by a Z-12 receiver.
%	      Typical call: z12head('blt03a97.149')

%Kai Borre 10-03-03
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2003/03/10 $

fidobs = fopen(observationfile);

% The size of a rawstructure header = 90 bytes
FA = fread(fidobs,42,'char');
FB = fread(fidobs,1,'int8');  
FC = fread(fidobs,1,'long');  % 4char
FD = fread(fidobs,43,'char');
rawheader = struct('version',setstr(FA(1:10)'), ...
    'raw_version',FA(11), ...
    'rcvr_type',setstr(FA(12:21)'), ...            
    'chan_ver',setstr(FA(22:31)'), ...
    'nav_ver',setstr(FA(32:41)'), ...
    'capability',FB(1), ...
    'reserved',FC(1), ... 
    'num_obs_types',FD(1), ...
    'spare',setstr(FD(2:43)'));

%%%%%%%% end z12head.m %%%%%%%%%%%%%%%%%%%%%%%%

