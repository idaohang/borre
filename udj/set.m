function sats(b);
% SATSUDJ’VNING
%
% Retningsobservationerne skrives som en matrix. Hver s›jle indeholder
% de m†lte v‘rdier for den enkelte retning i de forskellige satser.
% Nulretningen udelades.
% Svaret beskrives ved de udj‘vnede retninger, orienteringskonstanterne
% samt spredningen p† en retning m†lt med en sats.

% Skrevet af Kai Borre
% 3 januar 1994

    [s,r]=size(b);	     % s satser, r+1 retninger
    b=[zeros(1,s); b'];
    b=b(:);		     % b er observationsligningernes h›jreside
    r = r+1;
    A(1:r,1:r+s) = [eye(r) -ones(r,1) zeros(r,s-1)]; % s=1
    for i = 2:s-1, A((i-1)*r+1:i*r,1:r+s) = [
		   eye(r) zeros(r,i-1) -ones(r,1) zeros(r,s-i)]; end
    A((s-1)*r+1:r*s,1:r+s) = [eye(r) zeros(r,s-1) -ones(r,1)];
    B = A(:,1:r);
    C = A(:,r+1:r+s);
    BB = B*B'/s;
    CC = C*C'/r;
    E = eye(r*s);
    retninger = B'*b/s
    orienteringskonst = C'*(E-BB)*b/r
    sigma0 = sqrt(b'*(E-BB)*(E-CC)*b/((r-1)*(s-1)))
