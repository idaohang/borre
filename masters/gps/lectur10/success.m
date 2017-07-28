function P = success;
%SUCCESS Computation of success rate for given ambiguities

global D 
n = length(D);

P = 1;
for i = 1:n
    x = 1/(2*sqrt(D(i)));
    P = P*erf(x);
end
%%%%%%%%%%%%%%%% end success.m %%%%%%%%