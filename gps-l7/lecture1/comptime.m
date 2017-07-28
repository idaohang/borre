function comptime(bfile)
%COMPTIME Reads receiver clock offset from a binary Ashtech observation
%	       file and plots it.
%	       See also b_clock

%Kai Borre 03-22-97
%Copyright (c) by Kai Borre
%$Revision: 1.1 $  $Date: 1998/10/28  $

global outm

v_light = 299792458;
b_clock(bfile)
outm = outm(:,[1]);
Outm = outm/(v_light*1.e-9);
plot(Outm)
title('Receiver clock offset, a reset clock','Fontsize',16)
xlabel('Epochs, interval 15 s','Fontsize',16)
ylabel('Clock offset [ns]','Fontsize',16)
set(gca,'Fontsize',16)
print comptime -deps
%%%%%%%%%%%%%%%%% end comptime.m %%%%%%%%%%%%%%%%%%%
