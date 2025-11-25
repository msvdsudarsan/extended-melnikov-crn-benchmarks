function M = melnikov_substrate()
% Extended Melnikov functional for substrate inhibition CRN
params = struct("k", 0.82, "alpha", 1.50);

t = linspace(0,20,20000);
y = exp(-0.12.*t) .* cos(1.7.*t);
M = trapz(t, abs(y) .* params.k .* params.alpha);
end
