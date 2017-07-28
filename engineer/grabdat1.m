function Obs = grabdat1(fid, NoSv, NoObs)
%GRABDAT1  Positioned in a RINEX file at a selected epoch
%	       reads observations of NoSv satellites

%Kai Borre 03-24-02
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 2002/03/24  $

Obs = zeros(NoSv, NoObs);
if NoObs <= 5	  % This will typical be Turbo SII data
    for u = 1:NoSv
        lin = fgetl(fid);
        for k = 1:NoObs
            Obs(u,k) = str2num(lin(3+16*(k-1):16*k-2));
        end
    end
else		      % This will typical be Z12 or Javad data
    NoObs = 5;
    for u = 1:NoSv
        lin = fgetl(fid);
        lin_doppler = fgetl(fid);
        for k = 1:NoObs
            numm = str2num(lin(1+16*(k-1):16*k-2));
            if isempty(numm) == 1, Obs(u,k) = NaN; else Obs(u,k)= numm; end
        end
        for k = 1:2
            numd = str2num(lin_doppler(1+16*(k-1):16*k-2));
            if isempty(numd) == 1, Obs(u,5+k) = NaN; else Obs(u,5+k) = numd; end
        end    
    end
end
%%%%%%%%% end grabdat1.m %%%%%%%%%
