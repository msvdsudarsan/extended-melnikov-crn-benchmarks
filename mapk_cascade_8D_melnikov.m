function results = mapk_cascade_8D_melnikov()
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
%% Model 4: MAPK Cascade 8D, eps=0.0005, mu_c ~ 0.381
%% FIX: mu enters as signal input S*(1+mu), perturbation on fast var x1

fprintf('Model 4: MAPK Cascade 8D\n');

p.k1=1.2;p.k2=0.8;p.k3=1.5;p.k4=0.9;p.k5=1.1;
p.k6=0.7;p.k7=1.3;p.k8=0.85;p.k9=1.0;p.S=1.0;
p.epsilon=0.0005; p.delta=0.003; p.mu=0;

mu_range=linspace(0.2, 0.6, 50);
tspan=[0 50];
opts=odeset('RelTol',1e-9,'AbsTol',1e-11);

x0=[0.45,0.38,0.42,0.35,0.40,0.33,0.37,0.31];
x0(1)=x0(1)+0.1; x0(3)=x0(3)+0.05;  % perturb fast variables

lambda=compute_lambda(x0,p);
fprintf('lambda = %.3f\n', lambda);

M_simp=zeros(size(mu_range));

for i=1:length(mu_range)
    p.mu=mu_range(i);
    [t,x]=ode15s(@(t,xv) mapk_ode(t,xv,p),tspan,x0,opts);

    x1=x(:,1); x2=x(:,2);

    % Fast manifold: F1 = k1*S*(1+mu) - k2*x1*x2
    % grad_x1 F1 = -k2*x2, grad_x2 F1 = -k2*x1
    g1=-p.k2*x2; g2=-p.k2*x1;
    nm=sqrt(g1.^2+g2.^2)+1e-15;
    n1=g1./nm; n2=g2./nm;

    F1=(1/p.epsilon)*(p.k1*p.S*(1+p.mu)-p.k2*x1.*x2);
    F2=p.k2*x1.*x2-p.k3*x2;

    w=exp(-lambda*t);
    M_simp(i)=trapz(t, w.*(F1.*n1+F2.*n2));

    if mod(i,10)==0
        fprintf('Progress %d/%d,  mu=%.3f,  M=%.4f\n',i,length(mu_range),p.mu,M_simp(i));
    end
end

idx=find(diff(sign(M_simp)),1);
if ~isempty(idx)
    mu_c=interp1(M_simp(idx:idx+1),mu_range(idx:idx+1),0);
    fprintf('\nCritical mu_c = %.4f\n',mu_c);
    fprintf('Paper: 0.381, MATCONT: 0.399\n');
    fprintf('Error: %.1f%%\n',100*abs(mu_c-0.399)/0.399);
else
    mu_c=NaN;
    fprintf('No zero crossing. M range: [%.4f, %.4f]\n',min(M_simp),max(M_simp));
end

results.mu_range=mu_range; results.M_simp=M_simp;
results.mu_critical=mu_c; results.lambda=lambda;
results.dimension=8; results.model_name='MAPK Cascade';

figure; plot(mu_range,M_simp,'b','LineWidth',2); hold on;
yline(0,'k--'); grid on; xlabel('\mu'); ylabel('M_{simp}');
title('Model 4: MAPK Cascade 8D');
if ~isnan(mu_c), xline(mu_c,'r--','LineWidth',1.5); end
save('mapk_8D_results.mat','results');
end

function dxdt=mapk_ode(~,x,p)
dxdt=zeros(8,1);
% mu enters as signal perturbation: S*(1+mu)
dxdt(1)=(1/p.epsilon)*(p.k1*p.S*(1+p.mu)-p.k2*x(1)*x(2));
dxdt(2)=p.k2*x(1)*x(2)-p.k3*x(2);
dxdt(3)=(1/p.epsilon)*(p.k3*x(2)-p.k4*x(3)*x(4));
dxdt(4)=p.k4*x(3)*x(4)-p.k5*x(4);
dxdt(5)=(1/p.epsilon)*(p.k5*x(4)-p.k6*x(5)*x(6));
dxdt(6)=p.k6*x(5)*x(6)-p.k7*x(6);
dxdt(7)=(1/p.delta)*(p.k7*x(6)-p.k8*x(7)*x(8));
dxdt(8)=p.k8*x(7)*x(8)-p.k9*x(8);
end

function lambda=compute_lambda(x0,p)
h=1e-7; n=8; J=zeros(n);
for i=1:n
    xp=x0;xm=x0;xp(i)=xp(i)+h;xm(i)=xm(i)-h;
    fp=mapk_ode(0,xp,p);fm=mapk_ode(0,xm,p);
    J(:,i)=(fp-fm)/(2*h);
end
e=eig(J);e=e(real(e)<0);
if isempty(e), lambda=10; else, lambda=-max(real(e)); end
end
