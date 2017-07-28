xmax = 2;
steps = 200;
maxiter = 32;
Z = zeros(steps);

for m = 1:steps
    ci = (-xmax+2*xmax*m/steps);
    for n = 1:steps
        cr = -xmax+2*xmax*n/steps-.5;
        zr = cr;
        zi = ci;
        rmax = maxiter;
        for r = 0:maxiter
            zrn = zr*zr-zi*zi+cr;
            zin = 2*zi*zr+ci;
            zi = zin;
            zr = zrn;
            if (zr*zr+zi*zi) > 4, rmax = r;  break
            end
        end
        Z(m,n) = sqrt(rmax/maxiter);
    end
end