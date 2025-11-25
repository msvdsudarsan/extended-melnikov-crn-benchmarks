function dydt = oregonator_rhs(t,y)
x  = y(1);
y1 = y(2);
z  = y(3);

eps  = 0.05;
beta = 2;

dx  = (1/eps) * (y1 - x*y1 + x*(1 - x));
dy1 = x*y1 - beta*y1 + z;
dz  = eps * (x - z);

dydt = [dx; dy1; dz];
end
