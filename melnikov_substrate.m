function results = melnikov_substrate()
%% Paper Title: "An Efficient Adjoint-Free Melnikov-Lyapunov Diagnostic for
%%               Transition Detection in Slow-Fast Chemical Reaction Networks"
%% Author 1:    Sri Venkata Durga Sudarsan Madhyannapu
%% Author 2:    Pradheep Kumar S.
%%
%% Affiliation 1: Freshmen Engineering Department, NRI Institute of Technology
%%                (Autonomous), Pothavarappadu, Agiripalli, Vijayawada,
%%                521212, Andhra Pradesh, India
%% Affiliation 2: Research Scholar, Jawaharlal Nehru Technological University
%%                Kakinada, Andhra Pradesh, India
%% Affiliation 3: School of Basic Sciences, SRM University AP, Neerukonda,
%%                Mangalagiri Mandal, Guntur, 522240, Andhra Pradesh, India
%%
%% Journal:       Communications in Nonlinear Science and Numerical Simulation
%%                Elsevier, ISSN: 1007-5704
%% Manuscript ID: CNSNS-D-26-00848
%% Submitted:      19 February 2026
%% Status:        Under Review, 2026
%% SSRN:          https://ssrn.com/abstract=6275667
%% SSRN ID:        6275667 (Distributed: 02/20/2026)
%%
%% Paper: CNSNS Elsevier | Manuscript ID: CNSNS-D-26-00848
%% Model 1: Substrate Inhibition Oscillator (2D), eps=0.005
%% Algorithm: For each mu, integrate TRANSIENT from x0 near equilibrium.
%% M(mu) = integral w(t)*<F(phi_mu,mu), n(phi_mu)> dt
%% Sign of M changes at Hopf bifurcation mu_c.
%% Expected: mu_c ~ 0.438, MATCONT ref = 0.472 (paper Table 1)

fprintf('Model 1: Substrate Inhibition Oscillator (2D)\n');

p.epsilon = 0.005;
p.lambda  = 2.50;
tspan     = [0, 25];
mu_range  = linspace(0.10, 0.80, 60);

% Start NEAR but not AT equilibrium (x1_eq ~ 0.45, x2_eq = hill(x1_eq))
x1_0 = 0.55;  % perturbed from equilibrium
x0   = [x1_0; x1_0/(1+x1_0^2) + 0.05];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.02);

M_simp = zeros(size(mu_range));

for i = 1:length(mu_range)
    p.mu = mu_range(i);
    [t, x] = ode15s(@(t,x) sub_ode(t,x,p), tspan, x0, opts);

    x1 = x(:,1);  x2 = x(:,2);
    hill  = x1./(1+x1.^2);
    dhill = (1-x1.^2)./(1+x1.^2).^2;

    % Unit normal to slow manifold S: {x2 = x1/(1+x1^2)}
    % grad F = (-dhill, 1), normalise
    g1 = -dhill; g2 = ones(size(x1));
    nm = sqrt(g1.^2+g2.^2);
    n1 = g1./nm;  n2 = g2./nm;

    % Full RHS
    F1 = (1/p.epsilon)*(x2 - hill);
    F2 = 0.5 - x2 - hill + p.mu*x1;

    w = exp(-p.lambda*t);
    M_simp(i) = trapz(t, w.*(F1.*n1 + F2.*n2));

    if mod(i,10)==0
        fprintf('  Progress: %d/%d,  mu=%.3f,  M=%.4f\n', ...
            i,length(mu_range),mu_range(i),M_simp(i));
    end
end

idx = find(diff(sign(M_simp)),1);
if ~isempty(idx)
    mu_c = interp1(M_simp(idx:idx+1),mu_range(idx:idx+1),0);
    fprintf('\nDetected transition: mu_c = %.4f\n', mu_c);
    fprintf('MATCONT reference:   0.472  (paper Table 1)\n');
    fprintf('Relative error:      %.1f%%\n',100*abs(mu_c-0.472)/0.472);
else
    mu_c = NaN;
    fprintf('No zero crossing. M range: [%.4f, %.4f]\n',min(M_simp),max(M_simp));
    [~,imin]=min(abs(M_simp)); 
    fprintf('Closest to zero at mu=%.4f, M=%.6f\n',mu_range(imin),M_simp(imin));
end

results.mu_range=mu_range; results.M_simp=M_simp;
results.mu_critical=mu_c; results.lambda=p.lambda;
results.epsilon=p.epsilon; results.dimension=2;
results.model_name='Substrate Inhibition Oscillator 2D';

figure; plot(mu_range,M_simp,'b','LineWidth',2); hold on;
yline(0,'k--'); grid on; xlabel('\mu'); ylabel('M_{simp}');
title('Model 1: Substrate Inhibition 2D');
if ~isnan(mu_c), xline(mu_c,'r--','LineWidth',1.5); end
save('substrate_2D_results.mat','results');
end

function dxdt = sub_ode(~,x,p)
hill    = x(1)/(1+x(1)^2);
dxdt    = zeros(2,1);
dxdt(1) = (1/p.epsilon)*(x(2)-hill);
dxdt(2) = 0.5-x(2)-hill+p.mu*x(1);
end
