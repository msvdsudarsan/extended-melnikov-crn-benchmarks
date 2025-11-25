function dydt = substrate_rhs(t,y)
x = y(1); s = y(2);
dx = 2*x*(1 - x) - 0.45 * x * s;
ds = -0.3 * s + 0.1 * x^2;
dydt = [dx; ds];
end
