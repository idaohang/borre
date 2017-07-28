%REC_CLOC  Plotting receiver clock offsets
%	        from Turbo-SII and Z-12 receivers

%Kai Borre, March 18, 1997, at Columbus, Ohio
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date:1998/12/07  $

%Units are either seconds, meters, or radians
v_light = 299792458;     	% vacuum speed of light m/s
OffSet = kalclock('pta.96o','pta.nav',0);
fprintf('\nAdded %6.0f meters to obtain positive offsets\n',...
                                            -floor(min(OffSet)))
OffSet = OffSet - floor(min(OffSet));
t = 1:80;
z = 1:240;

%THE FOLLOWING IS DUE TO FIGURE PRODUCTION
%  We get a sample of Z-12 clock offset
%  The following line only to be run once---takes a minute or two
%  b_clock('b0810a94.076')

fid = fopen('clock.dat','r');
[T,count] = fread(fid,inf,'double');
fprintf('\nAdded %6.0f meters to obtain positive offsets\n\n',...
                        -floor(min(T(1:240))))
A = T(1:240)-floor(min(T(1:240)));
toc
figure;
semilogy(20*z, A/v_light,'b-',15*t,OffSet/v_light,'r--')
ylabel('Clock offset  [s]','Fontsize',16)
xlabel('Time  [s]','Fontsize',16)
set(gca,'FontSize',16)
hl = legend(' Clock with reset', ' Steered  clock ');
axes(hl);
set(gca,'FontSize',16)
refresh, pause
print rec_cloc -deps
fclose('all');
%%%%%%%%% end rec_cloc.m %%%%%%%%%
