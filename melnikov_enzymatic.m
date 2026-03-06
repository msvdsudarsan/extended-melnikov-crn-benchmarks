function results = melnikov_enzymatic()
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
%% Model 2: Enzymatic Feedback Network (3D), eps=0.008
%% Transient approach: M(mu) changes sign at Hopf bifurcation
%% Expected: mu_c ~ 0.326, MATCONT ref = 0.343

fprintf('Model 2: Enzymatic Feedback Network (3D)\n');

p.epsilon = 0.008;
p.lambda  = 3.20;
tspan     = [0, 30];
mu_range  = linspace(0.05, 0.70, 60);

% Start perturbed from equilibrium
y_eq = 1+sqrt(2);
x1_eq = 2*y_eq/(1+y_eq^2);
x0 = [x1_eq+0.2; y_eq+0.3; y_eq+0.2];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.05);

M_simp = zeros(size(mu_range));

for i = 1:length(mu_range)
    p.mu = mu_range(i);
    [t,x] = ode15s(@(t,x) enz_ode(t,x,p), tspan, x0, opts);

    x1=x(:,1); x2=x(:,2); x3=x(:,3);

    % Fast manifold: F_fast = -x1 + 2x2/(1+x2^2) = 0
    % grad in (x1,x2): g1=-1, g2=2(1-x2^2)/(1+x2^2)^2
    g1 = -ones(size(x1));
    g2 = 2*(1-x2.^2)./(1+x2.^2).^2;
    nm = sqrt(g1.^2+g2.^2)+1e-15;
    n1 = g1./nm;  n2 = g2./nm;

    F1 = (1/p.epsilon)*(-x1+2*x2./(1+x2.^2));
    F2 = x1 - x2./(1+x3) + p.mu*x2;

    w = exp(-p.lambda*t);
    M_simp(i) = trapz(t, w.*(F1.*n1 + F2.*n2));

    if mod(i,10)==0
        fprintf('  Progress: %d/%d,  mu=%.3f,  M=%.4f\n',...
            i,length(mu_range),mu_range(i),M_simp(i));
    end
end

idx = find(diff(sign(M_simp)),1);
if ~isempty(idx)
    mu_c = interp1(M_simp(idx:idx+1),mu_range(idx:idx+1),0);
    fprintf('\nDetected transition: mu_c = %.4f\n',mu_c);
    fprintf('MATCONT reference:   0.343\n');
    fprintf('Relative error:      %.1f%%\n',100*abs(mu_c-0.343)/0.343);
else
    mu_c = NaN;
    fprintf('No zero crossing. M range: [%.4f, %.4f]\n',min(M_simp),max(M_simp));
end

results.mu_range=mu_range; results.M_simp=M_simp;
results.mu_critical=mu_c; results.lambda=p.lambda;
results.epsilon=p.epsilon; results.dimension=3;
results.model_name='Enzymatic Feedback Network 3D';

figure; plot(mu_range,M_simp,'b','LineWidth',2); hold on;
yline(0,'k--'); grid on; xlabel('\mu'); ylabel('M_{simp}');
title('Model 2: Enzymatic Feedback 3D');
if ~isnan(mu_c), xline(mu_c,'r--','LineWidth',1.5); end
save('enzymatic_3D_results.mat','results');
end

function dxdt = enz_ode(~,x,p)
dxdt=zeros(3,1);
dxdt(1)=(1/p.epsilon)*(-x(1)+2*x(2)/(1+x(2)^2));
dxdt(2)=x(1)-x(2)/(1+x(3))+p.mu*x(2);
dxdt(3)=0.1*(x(2)-x(3));
end
