function results = brusselator_2D_melnikov()
%% brusselator_2D_melnikov.m
%% Model 7: Brusselator (2D), epsilon = 0.001
%% Reference: Prigogine & Lefever (1968)
%
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
%% Status:        Under Review, 2026

params.A = 1.0; params.B = 3.0; params.epsilon = 0.001;
mu_range = linspace(2.5, 3.5, 50);
M_simp   = zeros(size(mu_range));
tspan    = [0, 30];
opts     = odeset('RelTol',1e-9,'AbsTol',1e-11);
x0       = [params.A, params.B/params.A];
lambda   = compute_lambda_brus(params);
fprintf('lambda = %.3f\n', lambda);

for i = 1:length(mu_range)
    params.B = mu_range(i);
    [t,x]    = ode15s(@(t,x)brus_rhs(t,x,params), tspan, x0, opts);
    M_simp(i)= melnikov_integral(t, x, lambda);
    if mod(i,10)==0, fprintf('Progress %d/%d\n',i,length(mu_range)); end
end

idx = find(diff(sign(M_simp)), 1);
if isempty(idx), mu_c = NaN;
else, mu_c = interp1(M_simp(idx:idx+1), mu_range(idx:idx+1), 0); end

results.mu_range    = mu_range;
results.M_simp      = M_simp;
results.mu_critical = mu_c;
results.lambda      = lambda;
results.dimension   = 2;
results.model_name  = 'Brusselator';
fprintf('Critical mu_c = %.4f\n', mu_c);

figure;
plot(mu_range, M_simp, 'b', 'LineWidth', 2); hold on;
yline(0,'k--'); xlabel('B'); ylabel('M_{simp}');
title('Brusselator 2D – Adjoint-Free Melnikov Function'); grid on;
save('brusselator_2D_results.mat','results');
end

function dxdt = brus_rhs(~, x, p)
dxdt = zeros(2,1);
dxdt(1) = (1/p.epsilon)*(p.A - (p.B+1)*x(1) + x(1)^2*x(2));
dxdt(2) = p.B*x(1) - x(1)^2*x(2);
end

function lambda = compute_lambda_brus(p)
xeq = p.A; yeq = p.B/p.A;
J = [-(p.B+1)+2*xeq*yeq,  xeq^2;
      p.B - 2*xeq*yeq,   -xeq^2];
J(1,:) = J(1,:)/p.epsilon;
e = eig(J); e = e(real(e)<0);
if isempty(e), lambda = 3.5; else, lambda = -max(real(e)); end
end

function M = melnikov_integral(t, x, lambda)
p_vals = x(:,1);
dx = diff(x); dt = diff(t);
nrm = zeros(length(t)-1,1);
for i = 1:length(nrm), nrm(i) = norm(dx(i,:)/dt(i)); end
nrm = [nrm; nrm(end)];
w = exp(-lambda*(t - t(1)));
M = trapz(t, w.*p_vals.*nrm);
end
