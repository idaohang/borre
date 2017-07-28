%M12_DATA     Serial I/O from Motorola M12 GPS receiver
%             Download of date, time, position, speed,
%             satellite and channel data, oscillator,
%             UTC and GMT offsets, et cetera

%Kai Borre
%Copyright (c) by Kai Borre
%$Revision 1.0 $  $Date2001/10/31  $

% The actual setup uses the COM2 port
s = serial('COM3'); %%3
s.InputBufferSize = 50000;

%%for tt = 1:20
fopen(s);
fprintf(s,'*IDN?')
set(s,'BaudRate',9600);
set(s,'TimeOut',50); % 25
test = fread(s,4,'char');
if char(test)' == '@@Ha'
    
    %Date
    month = fread(s,1);
    day = fread(s,1);
    y = fread(s,2);
    % y = fread(s,1,'int16')
    year = y(1)*256+y(2);
    fprintf('\nDate: %2g/%2g/%4g',month,day,year) 
    
    %Time
    hour = fread(s,1,'uint8');
    minute = fread(s,1,'uint8');
    second = fread(s,1,'uint8');
    ns = fread(s,4,'uchar');
    nanosecond = ns(1)*256^3+ns(2)*256^2+ns(3)*256+ns(4);
    fprintf('\nTime: %2g:%2g:%11.9f\n',...
        hour,minute,second+nanosecond*10^(-9)) 
    
    %Position, i = 1 Filtered or Unfiltered Following Filter Select
    %Position, i = 2 Always Unfiltered
    for i = 1:2
        lat = fread(s,4,'uchar');
        lat_millisec = lat(1)*256^3+lat(2)*256^2+lat(3)*256+lat(4);
        latitude = rad2dms0(lat_millisec*pi/(180*60*60*1000));
        lon = fread(s,4,'uchar');
        lon_millisec = lon(1)*256^3+lon(2)*256^2+lon(3)*256+lon(4);
        longitude = rad2dms0(lon_millisec*pi/(180*60*60*1000));
        hei = fread(s,4,'uchar');
        height_GPS = (hei(1)*256^3+hei(2)*256^2+hei(3)*256+...
            hei(4))/100;
        hei = fread(s,4,'uchar');
        height_MSL = (hei(1)*256^3+hei(2)*256^2+hei(3)*256+...
            hei(4))/100;
        if i == 1
            fprintf('\nLatitude:     %02g %02g %05.3f',latitude)
            fprintf('\nLongitude: %3g %02g %05.3f',longitude)     
            fprintf('\nGPS Height:      %7.2f m\n',height_GPS)    
        end
    end
    
    %Speed/Heading
    VV = fread(s,2); % cm/s
    speed_3d = (VV(1)*256+VV(2))/100;
    fprintf('\n3D Speed:      %5.2f m/s',speed_3d)
    vv = fread(s,2); % cm/s
    speed_2d = (vv(1)*256+vv(2))/100;
    fprintf('\n2D Speed:      %5.2f m/s',speed_2d)
    hh = fread(s,2); % degree/10
    heading = (hh(1)*256+hh(2))/10;
    fprintf('\n2D Heading: %4.1f degrees\n',heading)
    
    %Geometry
    dd = fread(s,2); % dop/10
    fprintf('\nPDOP for 3D Fix, HDOP for 2D Fix, Zero for No Fix')
    DOP = (dd(1)*256+dd(2))/10;
    fprintf('\nCurrent DOP:   %3.1f',DOP)
    
    %Satellite Data
    PRN_vis = fread(s,1);
    fprintf('\nNumber of Visible Satellites:   %2g',PRN_vis)
    PRN_tracked = fread(s,1);
    fprintf('\nNumber of Tracked Satellites:  %2g\n',PRN_tracked)    
    
    %Channel Data
    for channel = 1:12   
        SVid(1,channel) = fread(s,1);
        mode(1,channel) = fread(s,1);
        signal_strength(1,channel) = fread(s,1); 
        IODE(1,channel) = fread(s,1);
        cs = fread(s,2);
        channel_status(1,channel) = cs(1)*256+cs(2);
    end
    fprintf('\n PRN       ')
    fprintf('%6g',SVid)    
    fprintf('\n Mode      ')
    fprintf('%6g',mode)
    fprintf('\n Strength')
    fprintf('%6g',signal_strength)
    fprintf('\n IODE    ')
    fprintf('%6g',IODE)
    fprintf('\n')
    
    ss = fread(s,2); %receiver status
    ssbit1 = bitget(ss(1),1:8);
    ssbit2 = bitget(ss(2),1:8); %not used
    pat = ssbit1(8:-1:6);
    pats = num2str(pat);
    switch pats
    case num2str([1 1 1]), strs = '3D Fix';    
    case num2str([1 1 0]), strs = '2D Fix';  
    case num2str([1 0 1]), strs = 'Propagate Mode';  
    case num2str([1 0 0]), strs = 'Position Hold';            
    case num2str([0 1 1]), strs = 'Acquiring Satellites';  
    case num2str([0 1 0]), strs = 'Bad Geometry';
    otherwise 
        strs = 'Unknown';    
    end
    fprintf('\nReceiver Status: %20s',strs)
    rr = fread(s,2); %reserved
    raw_range = rr(1)*256+rr(2);
    fprintf('\nRaw Range % 12.3f',raw_range); % Correct guess?
    
    %Oscillator and Clock Parameters
    cc = fread(s,2,'schar');
    c_bias = cc(1)*256+cc(2);
    fprintf('\nClock Bias:         %-5.0f ns',c_bias)
    oooo = fread(s,4);
    osc_offset = oooo(1)*256^3+oooo(2)*256^2+oooo(3)*256+oooo(4);
    fprintf('\nOscillator Offset: %6g Hz',osc_offset)
    TT = fread(s,2,'schar'); %half-degrees Celcius
    temperature = (TT(1)*256+TT(2))/2;
    fprintf('\nTemperature: %-2g C\n',temperature)
    
    %UTC parameters
    u = fread(s,1);
    UTC_o = bitget(u,1:6);
    UTC_offset = UTC_o(1)+UTC_o(2)*2+UTC_o(3)*2^2+UTC_o(4)*2^3+...
        UTC_o(5)*2^4;
    fprintf('\nUTC offset from GPS Time:  %+2g s',UTC_offset)
    
    %GMT Offset
    s0 = fread(s,1);
    GSM_offh = fread(s,1);
    if s0 == 'ff', GSM_offh = -GSM_offh; end  % Check if correct
    GSM_offm = fread(s,1);
    fprintf('\nGMT Offset:     %+1g:%02g h\n',GSM_offh,GSM_offm)
    vvvvvv = fread(s,6,'char');
    fprintf('\nID Tag: %6s\n\n',char(vvvvvv'))
    check_sum = fread(s,1);
end
fclose(s)
%%end
delete(s)
clear s
%%%%%%%%%%%%%%% end m12_data.m  %%%%%%%%%%%%%%

