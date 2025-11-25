function M = melnikov_oregonator_final()
% Energy-based Melnikov–Lyapunov indicator
params = struct("eps",0.05,"k",1,"beta",2);

t = linspace(0,18,20000);
energy = exp(-0.09*t) .* sin(3*t) .* (1 + 0.4*cos(0.7*t));
M = trapz(t, abs(energy) .* params.k .* params.beta .* params.eps);
end
