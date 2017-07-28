xmax = 2;
steps = 200;
maxiter = 32;
Z = 0;
i = sqrt(-1);

for m = 1:steps
    for n = 1:steps
        c = -xmax+2*xmax*n/steps-.5+i*(-xmax+2*xmax*m/steps);
        z = c;
        for r = 0:maxiter
            z = z*z+c;
            if abs(z) > 2, break
            end
        end
        Z(m,n) = sqrt(r/maxiter);
    end
end