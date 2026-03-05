function results = glycolysis_7D_melnikov()
%% glycolysis_7D_melnikov.m
%% Model 3: Extended Glycolysis (7D), epsilon = 0.0008
%% Reference: Goldbeter (1996) Biochemical Oscillations and Cellular Rhythms
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

params.v0=0.36; params.k1=1.0; params.k2=6.0; params.k3=0.8;
params.k4=1.2;  params.k5=0.9; params.k6=1.5; params.k7=0.7;
params.k8=1.1;  params.epsilon=0.0008; params.delta=0.005;

mu_range = linspace(0.1, 0.5, 50);
M_simp   = zeros(size(mu_range));
tspan    = [0, 40];
options  = odeset('RelTol',1e-9,'AbsTol',1e-11);
x0       = [0.35, 0.28, 0.42, 0.31, 0.38, 0.29, 0.33];
lambda   = compute_lambda_glycolysis(x0, params);
fprintf('lambda = %.2f\n', lambda);

for i = 1:length(mu_range)
    params.mu = mu_range(i);
    [t, x]    = ode15s(@(t,x)glycolysis_rhs(t,x,params), tspan, x0, options);
    M_simp(i) = melnikov_integral_glycolysis(t, x, lambda);
    if mod(i,10)==0, fprintf('Progress: %d/%d\n',i,length(mu_range)); end
end

zi = find(diff(sign(M_simp)), 1);
if ~isempty(zi)
    mu_c = interp1(M_simp(zi:zi+1), mu_range(zi:zi+1), 0);
else, mu_c = NaN; end

results.mu_range    = mu_range;
results.M_simp      = M_simp;
results.mu_critical = mu_c;
results.lambda      = lambda;
results.dimension   = 7;
results.model_name  = 'Glycolysis 7D';
fprintf('Critical mu_c = %.4f\n', mu_c);

figure;
plot(mu_range, M_simp, 'b', 'LineWidth', 2); hold on;
yline(0,'k--'); xlabel('\mu'); ylabel('M_{simp}');
title('Glycolysis 7D – Adjoint-Free Melnikov Function'); grid on;
save('glycolysis_7D_results.mat','results');
end

function dxdt = glycolysis_rhs(~, x, p)
dxdt    = zeros(7,1);
dxdt(1) = (1/p.epsilon)*(p.v0 - (p.k1*x(1)*x(2)^2)/(1+x(2)^2));
dxdt(2) = (p.k1*x(1)*x(2)^2)/(1+x(2)^2) - p.k2*x(2) - p.k3*x(2)*x(3);
dxdt(3) = p.k3*x(2)*x(3) - p.k4*x(3);
dxdt(4) = (1/p.delta)*(p.k4*x(3) - p.k5*x(4)*x(5));
dxdt(5) = p.k5*x(4)*x(5) - p.k6*x(5);
dxdt(6) = p.k6*x(5) - p.k7*x(6)*x(7);
dxdt(7) = (1/p.delta)*(p.k7*x(6)*x(7) - p.k8*x(7));
end

function lambda = compute_lambda_glycolysis(x0, p)
h = 1e-6; J = zeros(7);
for i = 1:7
    xp = x0; xm = x0; xp(i)=xp(i)+h; xm(i)=xm(i)-h;
    fp = glycolysis_rhs(0,xp,p); fm = glycolysis_rhs(0,xm,p);
    J(:,i) = (fp-fm)/(2*h);
end
e = eig(J); e = e(real(e)<0);
if isempty(e), lambda = 5; else, lambda = -max(real(e)); end
end

function M = melnikov_integral_glycolysis(t, x, lambda)
p_vals = x(:,2);
dx = diff(x(:,1:3)); dt = diff(t);
nrm = zeros(length(t)-1,1);
for i = 1:length(nrm), nrm(i) = norm(dx(i,:)/dt(i)); end
nrm = [nrm; nrm(end)];
w = exp(-lambda*(t - t(1)));
M = trapz(t, w.*p_vals.*nrm);
end
