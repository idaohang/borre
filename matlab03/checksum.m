function C = checksum(let1,let2,par)
%CHECKSUM  Generating Checksum for M12 Motorola GPS card
%          C is an integer

%Written By Kai Borre 
%June 3, 2001

a = dec2bin(let1,8);
b = dec2bin(let2,8);
xxr = a ~= b;
c = dec2bin(par,8);
xxxr = mod(xxr+c,2);
result = sprintf('%1g',xxxr);
C = bin2dec(result);
%char(cp2) % ASCII character
%%%%%%%%%%%%%%%%% end checksum.m  %%%%%%%%%%%%%%%5
