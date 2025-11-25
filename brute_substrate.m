function R = brute_substrate()
% High-accuracy reference simulation for error metric
opts = odeset("RelTol",1e-9,"AbsTol",1e-10);
[~,X] = ode15s(@substrate_rhs,[0 20],[0.9 0.1],opts);
R = max(abs(X(:,1)));
end
