%Partielle differentialkvotienter af afstandsobservationen
echo on
syms s Yj Yk Xj Xk real;
b = sqrt((Yk-Yj)^2+(Xk-Xj)^2);
J = jacobian(b,[Yj Yk Xj Xk]);
J = subexpr(J,'s2')
J = simplify(J)
J = collect(J)
echo off
%%%%%%%%%%%%%%%%% end part.m  %%%%%%%%%%%%%%%%%