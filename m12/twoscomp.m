function r = twoscomp(a)
%TWOSCOMP  Two's complement
%          In DSP, variables are often represented as fixed-point,
%          two's complement fractions. In this representation, the 
%          binary point is to the right of the MSB (most significant bit) 
%          which is also the sign bit (1 is minus, 0 is plus).
%          Two's complement positive numbers are in the natural binary
%          form. A negative number is formed from the corresponding 
%          positive number by complementing all the bits of the positive 
%          number and then adding 1 LSB (least significant bit).

%Kai Borre
%Copyright (c) by Kai Borre
%$Revision 1.0 $  $Date2001/10/31  $

col = size(a,2);
if char(a(1)) == '1'   % minus sign
    v = a-'0';         % Convert to numbers
    I = find(v == 1);
    v = ones(1,col);
    v(I) = 0;
    vec = dec2bin(v)';
    r = -(bin2dec(vec)+1);
else
    r = bin2dec(a);    
end
%%%%%%%%%%%%%%%%%%%% end twoscomp.m  %%%%%%%%%%%%%%%%%%