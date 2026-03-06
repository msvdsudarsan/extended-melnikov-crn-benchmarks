function results = glycolysis_7D_melnikov()
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
%% Model 3: Glycolysis 7D, eps=0.0008, mu_c ~ 0.267
%% FIX: Perturbation enters x1 (fast substrate), not x2 or x5
%% dF/dmu = (x1, 0, 0, 0, 0, 0, 0)^T when mu*x1 added to x1 eq
%% Alternatively: mu enters as input flux v0 -> v0*(1+mu)

fprintf('Model 3: Glycolysis 7D\n');

p.v0=0.36; p.k1=1.0; p.k2=6.0; p.k3=0.8;
p.k4=1.2;  p.k5=0.9; p.k6=1.5; p.k7=0.7; p.k8=1.1;
p.epsilon=0.0008; p.delta=0.005; p.mu=0;

mu_range = linspace(0.10, 0.50, 50);
tspan    = [0, 40];
opts     = odeset('RelTol',1e-9,'AbsTol',1e-11);

% Start perturbed from approximate equilibrium
x0 = [0.35,0.28,0.42,0.31,0.38,0.29,0.33];
x0(1) = x0(1)+0.1; x0(2) = x0(2)+0.05;  % perturb fast variables

lambda = compute_lambda(x0, p);
fprintf('lambda = %.2f\n', lambda);

M_simp = zeros(size(mu_range));

for i = 1:length(mu_range)
    p.mu = mu_range(i);
    [t, x] = ode15s(@(t,xv) gly_ode(t,xv,p), tspan, x0, opts);

    x1=x(:,1); x2=x(:,2);

    % Normal to fast manifold: F_fast = v0*(1+mu) - k1*x1*x2^2/(1+x2^2)
    % grad_x1 F = -k1*x2^2/(1+x2^2)
    % grad_x2 F = -k1*x1*2x2/(1+x2^2)^2
    g1 = -p.k1*x2.^2./(1+x2.^2);
    g2 = -p.k1*x1.*2.*x2./(1+x2.^2).^2;
    nm = sqrt(g1.^2+g2.^2)+1e-15;
    n1 = g1./nm;  n2 = g2./nm;

    F1 = (1/p.epsilon)*(p.v0*(1+p.mu)-(p.k1*x1.*x2.^2)./(1+x2.^2));
    F2 = (p.k1*x1.*x2.^2)./(1+x2.^2)-p.k2*x2-p.k3*x2.*x(:,3);

    w = exp(-lambda*t);
    M_simp(i) = trapz(t, w.*(F1.*n1 + F2.*n2));

    if mod(i,10)==0
        fprintf('Progress: %d/%d,  mu=%.3f,  M=%.4f\n',i,length(mu_range),p.mu,M_simp(i));
    end
end

idx = find(diff(sign(M_simp)),1);
if ~isempty(idx)
    mu_c = interp1(M_simp(idx:idx+1),mu_range(idx:idx+1),0);
    fprintf('\nCritical mu_c = %.4f\n',mu_c);
    fprintf('Paper: 0.267, MATCONT: 0.283\n');
    fprintf('Error: %.1f%%\n',100*abs(mu_c-0.283)/0.283);
else
    mu_c = NaN;
    fprintf('No zero crossing. M range: [%.4f, %.4f]\n',min(M_simp),max(M_simp));
end

results.mu_range=mu_range; results.M_simp=M_simp;
results.mu_critical=mu_c; results.lambda=lambda;
results.dimension=7; results.model_name='Glycolysis 7D';

figure; plot(mu_range,M_simp,'b','LineWidth',2); hold on;
yline(0,'k--'); grid on; xlabel('\mu'); ylabel('M_{simp}');
title('Model 3: Glycolysis 7D');
if ~isnan(mu_c), xline(mu_c,'r--','LineWidth',1.5); end
save('glycolysis_7D_results.mat','results');
end

function dxdt = gly_ode(~,x,p)
dxdt=zeros(7,1);
% mu enters as input flux perturbation: v0*(1+mu)
dxdt(1)=(1/p.epsilon)*(p.v0*(1+p.mu)-(p.k1*x(1)*x(2)^2)/(1+x(2)^2));
dxdt(2)=(p.k1*x(1)*x(2)^2)/(1+x(2)^2)-p.k2*x(2)-p.k3*x(2)*x(3);
dxdt(3)=p.k3*x(2)*x(3)-p.k4*x(3);
dxdt(4)=(1/p.delta)*(p.k4*x(3)-p.k5*x(4)*x(5));
dxdt(5)=p.k5*x(4)*x(5)-p.k6*x(5);
dxdt(6)=p.k6*x(5)-p.k7*x(6)*x(7);
dxdt(7)=(1/p.delta)*(p.k7*x(6)*x(7)-p.k8*x(7));
end

function lambda=compute_lambda(x0,p)
h=1e-6; J=zeros(7);
for i=1:7
    xp=x0; xm=x0; xp(i)=xp(i)+h; xm(i)=xm(i)-h;
    fp=gly_ode(0,xp,p); fm=gly_ode(0,xm,p);
    J(:,i)=(fp-fm)/(2*h);
end
e=eig(J); e=e(real(e)<0);
if isempty(e), lambda=5; else, lambda=-max(real(e)); end
end
