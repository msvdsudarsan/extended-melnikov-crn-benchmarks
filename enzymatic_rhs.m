function dydt = enzymatic_rhs(t,y)
x  = y(1);
y1 = y(2);
z  = y(3);

dx  = 1.5*x*(1 - x) - 0.4*x*z;
dy1 = 0.2*x - 0.8*y1;
dz  = 0.5*y1 - 0.6*z;

dydt = [dx; dy1; dz];
end
