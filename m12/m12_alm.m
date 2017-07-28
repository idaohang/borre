%M12_ALM    Serial I/O from Motorola M12 GPS board.
%           Downloading an almanac by sending the @@Be option.
%           Data are formatted according to the document:
%              Global Positioning System, Standard Positioning
%              Service, Signal Specification, 2nd Edition, June 2, 1995.
%              Section 2.4

%Kai Borre
%Copyright (c) by Kai Borre
%$Revision 1.0 $  $Date2001/10/31  $

Pi = 3.1415926535898; % pi value for GPS
s = serial('COM3');
s.OutputBufferSize = 512;
s.InputBuffersize = 50000;
fopen(s);
s.BaudRate = 9600;
set(s,'TimeOut',1);
s.RecordMode = 'index';
s.RecordDetail = 'verbose';
s.RecordName = 'motor.tex';
fprintf(s,'*IDN?')
C = checksum('B','e',0);
%C = checksum('G','j',0);
% Sending request for almanac data
fwrite(s,['@@Be' char(0) char(C) char(13) char(10)],'uint8') 
%fwrite(s,['@@Gj' char(C) char(13) char(10)],'uint8') 

t = 0;
while 1
    at = fread(s,2,'uint8');
    if char(at)' == '@@'
        letters = fread(s,2,'uint8');
        if char(letters)' == 'Cb'     
            t = t+1;
            sp = fread(s,2,'uint8');    
            subframe = sp(1);
            page = sp(2);
            fprintf('\n\n Subframe %2g Page %2g\n',subframe,page)
            for i = 3:10
                word(1:3,i) = fread(s,3,'uint8'); 
            end
            % Unpacking and Assigning the Kepler Elements
            svb = dec2bin(word(1,3),8);
            nosv = bin2dec(svb(1:2));
            sv = bin2dec(svb(3:8));
            if nosv ~= 1
                disp('No valid almanac')
            end
            
            % Eccentricity, e
            e = dec2bin(word(2:3,3),8); 
            ecc = bin2dec([e(1,:),e(2,:)])*2^(-21);
           
            % Time of almanac, t_0a
            to = dec2bin(word(1,4));
            toa = bin2dec(to)*2^12;            
         
            % Inclination, i
            d_ib = dec2bin(word(2:3,4),8);  
            delta_i = twoscomp([d_ib(1,:),d_ib(2,:)])*2^(-19);            
            i = (0.3+delta_i)*Pi;
            
            % Rate of right ascension, Omegadot
            Odot_b = dec2bin(word(1:2,5),8);   
            Omegadot = twoscomp([Odot_b(1,:),Odot_b(2,:)])*2^(-38)*Pi;           
                 
            % Health
            h_b = dec2bin(word(3,5),8);   
            health = h_b(1:3); %the 3 most significant bits, cf. Table 2-9
           
            % Square root of semi-major axis, sqrt(a)
            % Table 2-8 is erroneous. # of bits shall be 24, and not 16
            s_Ab = dec2bin(word(1:3,6),8);
            sqrt_A = bin2dec([s_Ab(1,:),s_Ab(2,:),s_Ab(3,:)])*2^(-11);           
               
            % Right Ascension, Omega_0
            O_0b = dec2bin(word(1:3,7),8);
            Omega_0 = twoscomp([O_0b(1,:),O_0b(2,:),O_0b(3,:)])*2^(-23)*Pi;
         
            % Argument of Perigee, omega
            o_b = dec2bin(word(1:3,8),8);
            omega = twoscomp([o_b(1,:),o_b(2,:),o_b(3,:)])*2^(-23)*Pi;
        
            % Mean Anomaly, mu_0
            m_0b = dec2bin(word(1:3,9),8);
            mu_0 = twoscomp([m_0b(1,:),m_0b(2,:),m_0b(3,:)])*2^(-23)*Pi;
        
            % Satellite clock offset, af_0
            af_0b = dec2bin(word(1,10),8);
            part = dec2bin(word(3,10),8);
            af_0bb = part(4:6);
            af_0 = twoscomp([af_0b,af_0bb])*2^(-20);
            
            % Rate of satellite clock offset, af_1
            af_1b = dec2bin(word(2,10),8);
            af_1bb = dec2bin(word(3,10),8); 
            af_1 = twoscomp([af_1b,af_1bb(1:3)])*2^(-38);
     
            % Week Number, WN
            if t == 33
                sp = dec2bin(word(1,3),8);   
                subframe = bin2dec(sp(1:2));
                page = bin2dec(sp(3:8));
                fprintf('\n\n Subframe %2g Page %2g\n',subframe,page)
                t_oab = dec2bin(word(2,3));
                t_0 = bin2dec(t_oab);
                wn_b = dec2bin(word(3,3)); 
                wn = bin2dec(wn_b);
                fprintf('\n Week Number:   %4.0f\n',wn)    
            end 
            checksum = fread(s,1,'uint8');
            carriage = fread(s,1,'uint8');
            linefeed = fread(s,1,'uint8');
            if t == 34, break, end
            
            fprintf('\n ID:                                \t\t\t %2g',sv)
            fprintf('\n Health (3 MSB):                    \t\t\t\t\t %3s',health)
            fprintf('\n Eccentricity                       \t %0.10e',ecc)
            fprintf('\n Time of Applicability (sec):       \t %8.0f',toa)
            fprintf('\n Orbital Inclination (rad):         \t %0.10e',i)
            fprintf('\n Rate of Right Ascen (rad/s):       \t %0.10e',Omegadot)     
            fprintf('\n Sqrt of Semi-Major Axis (m^(1/2)): \t %10.6f',sqrt_A)
            fprintf('\n Right Ascension (rad):             \t %0.10e',Omega_0)
            fprintf('\n Argument of Perigee (rad):         \t %0.10e',omega)
            fprintf('\n Mean Anomaly (rad):                \t %0.10e',mu_0)
            fprintf('\n Satellite Clock Offset (s):        \t %0.10e',af_0)
            fprintf('\n Rate of Satellite Clock Offset (s/s): \t %0.10e',af_1)
        end 
    end
end        
fclose(s)
delete(s)
clear s
%%%%%%%%%%%%% end m12_alm.m  %%%%%%%%%%%%%%

