function tR = brute_oregonator_final()
% Reference runtime only (state-based error not meaningful)
opts = odeset("RelTol",1e-9,"AbsTol",1e-10);
tic;
ode15s(@oregonator_rhs,[0 2],[1 1 1],opts);
tR = toc;
end
