function corrcols = repair_c(prc,prcc);

spikes = diff(prc);
j = find(abs(spikes) > 280000);
for k = 1:size(j)
   if spikes(j) < 0
      corr = 299792.458;
   else
      corr = -299792.458;
   end
   for l = j(k)+1:size(prc,1)
      prc(l,1) = prc(l,1)+corr;
      if nargin > 1
	 prcc(l,1) = prcc(l,1)+corr;
      end
  end
end
corrcols = prc;
if nargin > 1
   corrcols = [prc prcc];
end
%%%%%%%%%%%%%%%%%%%%%%%% end repair_c.m  %%%%%%%%%%%%%%%%%
