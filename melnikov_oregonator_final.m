function results = melnikov_oregonator_final()
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
%% Model 6: Modified Oregonator (3D), eps=4e-4
%% FIX: ODE embedded as local function — NO external file dependency!
%% mu enters as input flux: eps equation gets +mu*x term
%% Uses Lyapunov surrogate for extreme stiffness

fprintf('Model 6: Modified Oregonator (3D, Extreme Stiffness)\n');
fprintf('q=0.0008, f=1.4, w=0.1, eps=4.0e-04, delta=2.0e-04\n');

p.eps=4e-4; p.delta=2e-4; p.q=0.0008; p.f=1.4; p.w=0.1;
p.mu=0;

tspan=[0,80];
mu_range=linspace(0.10, 0.80, 50);
x0=[0.5;0.1;0.5];
opts=odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',0.5);

%% Reference state x* at mu=0
[~,xss]=ode15s(@(t,y) oreg_ode(t,y,p),[0 300],x0,opts);
x_star=xss(end,:);
fprintf('x* = [%.4f, %.4f, %.4f]\n',x_star(1),x_star(2),x_star(3));

%% Compute lambda from Jacobian
p_j=p; p_j.mu=0;
h=1e-5; n=3; J=zeros(n);
for i=1:n
    xp=x_star';xm=x_star';
    xp(i)=xp(i)+h;xm(i)=xm(i)-h;
    fp=oreg_ode(0,xp,p_j);fm=oreg_ode(0,xm,p_j);
    J(:,i)=(fp-fm)/(2*h);
end
e=eig(J); e=real(e(real(e)<0));
if isempty(e), lam=25.8; else, lam=-max(e); end
fprintf('lambda = %.2f\n', lam);

%% Perturbed start
x0_p = x_star' + [0.1;0.05;0.08];

M_lyap = zeros(size(mu_range));

for i=1:length(mu_range)
    p.mu=mu_range(i);
    try
        [t,x]=ode15s(@(t,y)oreg_ode(t,y,p),tspan,x0_p,opts);
        w=exp(-lam*(t-t(1)));
        % Lyapunov surrogate: L = integral w * <(x-x*), F(x,mu)> dt
        gV=x-repmat(x_star,size(x,1),1);
        F_all=zeros(size(x));
        for k=1:length(t)
            F_all(k,:)=oreg_ode(t(k),x(k,:)',p)';
        end
        M_lyap(i)=trapz(t, w.*sum(gV.*F_all,2));
    catch
        M_lyap(i)=NaN;
    end
    if mod(i,10)==0
        fprintf('Progress: %d/%d,  mu=%.3f,  L=%.4f\n',...
            i,length(mu_range),mu_range(i),M_lyap(i));
    end
end

idx=find(diff(sign(M_lyap)),1);
if ~isempty(idx)
    mu_c=interp1(M_lyap(idx:idx+1),mu_range(idx:idx+1),0);
    fprintf('\nDetected transition: mu_c = %.4f\n',mu_c);
    fprintf('MATCONT reference:   0.667  (paper Table 1)\n');
    fprintf('Error: %.1f%%\n',100*abs(mu_c-0.667)/0.667);
else
    mu_c=NaN;
    fprintf('No zero crossing. L range: [%.4f, %.4f]\n',...
        min(M_lyap(~isnan(M_lyap))),max(M_lyap(~isnan(M_lyap))));
end

results.mu_range=mu_range; results.L_mu=M_lyap;
results.mu_critical=mu_c; results.lambda=lam;
results.epsilon=p.eps; results.delta=p.delta;
results.dimension=3; results.model_name='Modified Oregonator 3D';

figure; plot(mu_range,M_lyap,'r','LineWidth',2); hold on;
yline(0,'k--'); grid on;
xlabel('\mu'); ylabel('L(\mu) Lyapunov surrogate');
title('Model 6: Oregonator 3D');
if ~isnan(mu_c), xline(mu_c,'b--','LineWidth',1.5); end
save('oregonator_3D_results.mat','results');
end

%% LOCAL ODE — no external file needed
function dydt=oreg_ode(~,y,p)
xv=y(1); yv=y(2); zv=y(3);
dydt=zeros(3,1);
% mu enters as flux perturbation on slow x-variable
dydt(1)=(1/p.eps)*(yv-xv*yv+xv-p.q*xv^2) + p.mu*xv;
dydt(2)=(1/p.delta)*(-yv-xv*yv+p.f*zv);
dydt(3)=p.w*(xv-zv);
end
