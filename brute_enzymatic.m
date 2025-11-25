function R = brute_enzymatic()
opts = odeset("RelTol",1e-9,"AbsTol",1e-10);
[~,X] = ode15s(@enzymatic_rhs,[0 30],[1.0 0.2 0.3],opts);
R = max(abs(X(:,1)));
end
