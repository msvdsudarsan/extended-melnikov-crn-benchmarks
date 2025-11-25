function M = melnikov_enzymatic()
% Extended Melnikov functional for enzymatic feedback CRN
params = struct("k", 1.2, "beta", 0.85, "gamma", 0.60);

t = linspace(0,30,20000);
weight = exp(-0.08*t);
core = sin(1.15*t) + 0.4*cos(0.65*t);
M = trapz(t, weight .* abs(core) .* params.k .* params.beta .* params.gamma);
end
