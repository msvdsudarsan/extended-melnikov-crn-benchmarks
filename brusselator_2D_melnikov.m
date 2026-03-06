function results = brusselator_2D_melnikov()
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
%% Model 7: Brusselator (2D), A=1, eps=0.001
%% B is bifurcation parameter. Hopf at B_c ~ 2.917.
%% Expected: mu_critical (=B_c) ~ 2.917, MATCONT ref = 3.000

p.A = 1.0;  p.epsilon = 0.001;

mu_range = linspace(2.0, 4.0, 60);
M_simp   = zeros(size(mu_range));

% Integration window — transient from perturbed start
tspan = [0, 30];
opts  = odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',0.05);

lambda = 3.5;
fprintf('lambda = %.3f\n', lambda);

for i = 1:length(mu_range)
    p.B = mu_range(i);
    x_eq = p.A;  y_eq = p.B/p.A;
    x0 = [x_eq + 0.3; y_eq + 0.3];  % perturbed from equilibrium

    [t, x] = ode15s(@(t,xv) brus_ode(t,xv,p), tspan, x0, opts);

    x1=x(:,1);  x2=x(:,2);

    % Fast manifold: A-(B+1)*x1+x1^2*x2=0
    % grad: g1=-(B+1)+2x1x2, g2=x1^2
    g1 = -(p.B+1)+2*x1.*x2;
    g2 = x1.^2;
    nm = sqrt(g1.^2+g2.^2)+1e-15;
    n1 = g1./nm;  n2 = g2./nm;

    F1 = (1/p.epsilon)*(p.A-(p.B+1)*x1+x1.^2.*x2);
    F2 = p.B*x1 - x1.^2.*x2;

    w = exp(-lambda*t);
    M_simp(i) = trapz(t, w.*(F1.*n1 + F2.*n2));

    if mod(i,10)==0
        fprintf('Progress %d/%d,  B=%.2f,  M=%.4f\n',i,length(mu_range),p.B,M_simp(i));
    end
end

idx = find(diff(sign(M_simp)),1);
if ~isempty(idx)
    mu_c = interp1(M_simp(idx:idx+1),mu_range(idx:idx+1),0);
    fprintf('\nCritical B_c = %.4f\n',mu_c);
    fprintf('Paper: 2.917, MATCONT ref: 3.000\n');
    fprintf('Error vs MATCONT: %.1f%%\n',100*abs(mu_c-3.000)/3.000);
else
    mu_c = NaN;
    fprintf('No zero crossing. M range: [%.4f, %.4f]\n',min(M_simp),max(M_simp));
end

results.mu_range=mu_range; results.M_simp=M_simp;
results.mu_critical=mu_c; results.lambda=lambda;
results.dimension=2; results.model_name='Brusselator';

figure; plot(mu_range,M_simp,'b','LineWidth',2); hold on;
yline(0,'k--'); grid on;
xlabel('B'); ylabel('M_{simp}');
title('Brusselator 2D Melnikov Function');
if ~isnan(mu_c), xline(mu_c,'r--','LineWidth',1.5); end
save('brusselator_2D_results.mat','results');
end

function dxdt = brus_ode(~,x,p)
dxdt=zeros(2,1);
dxdt(1)=(1/p.epsilon)*(p.A-(p.B+1)*x(1)+x(1)^2*x(2));
dxdt(2)=p.B*x(1)-x(1)^2*x(2);
end
