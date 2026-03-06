function results = circadian_10D_melnikov()
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
%% Model 5: Circadian Rhythm 10D, eps=0.0002
%% FIX 1: lambda from Jacobian (NOT 1/eps — that kills the integral!)
%% FIX 2: mu enters as vs*(1+mu) — input flux perturbation on fast x1
%% Uses Lyapunov surrogate for extreme stiffness
%% Expected: mu_c ~ 0.523, MATCONT ref = 0.556

fprintf('Model 5: Circadian Rhythm 10D\n');

p.vs=1.6;p.vm=0.8;p.Km=0.5;p.ks=0.6;p.vd=0.95;p.Kd=0.3;
p.k1=0.4;p.k2=0.2;p.KI=1.0;p.n=4;p.v1=0.7;p.K1=0.5;
p.k3=0.4;p.k4=0.2;
p.vs_p=1.5;p.vm_p=0.75;p.Km_p=0.5;p.ks_p=0.65;p.vd_p=1.0;p.Kd_p=0.35;
p.k1_p=0.45;p.k2_p=0.22;p.v1_p=0.75;p.K1_p=0.55;p.k3_p=0.42;p.k4_p=0.21;
p.epsilon=0.0002; p.delta=0.002; p.mu=0;

mu_range=linspace(0.3,0.7,40);
tspan=[0 60];
opts=odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',0.1);
x0=[0.5,0.4,0.45,0.35,0.42,0.48,0.38,0.44,0.36,0.40];
x0_p=x0; x0_p(1)=x0_p(1)+0.1; x0_p(6)=x0_p(6)+0.1;  % perturb

%% Compute lambda from Jacobian (NOT 1/epsilon!)
lambda = compute_lambda_circ(x0, p);
fprintf('lambda = %.2f  (from Jacobian, not 1/eps)\n', lambda);

%% Reference state x* for Lyapunov surrogate
lopts=odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',0.1);
[~,xss]=ode15s(@(t,x)circ_ode(t,x,p),[0 300],x0,lopts);
x_star=xss(end,:);

M_simp=zeros(size(mu_range));
M_lyap=zeros(size(mu_range));

for i=1:length(mu_range)
    p.mu=mu_range(i);
    try
        [t,x]=ode15s(@(t,xv)circ_ode(t,xv,p),tspan,x0_p,opts);
        w=exp(-lambda*(t-t(1)));

        % M_simp: fast manifold projection (x1,x5,x10 are fast)
        % Fast eq x1: F1=vs*(1+mu)*KI^n/(KI^n+x9^n)-vm*x1/(Km+x1)
        % grad_x1 F1 = -vm/(Km+x1)^2; grad_x9 F1 = -n*vs...
        x1=x(:,1); x9=x(:,9);
        Kn=p.KI^p.n;
        gF1 = -p.vm./(p.Km+x1).^2;
        gF9 = -p.n*p.vs*(1+p.mu)*Kn.*x9.^(p.n-1)./(Kn+x9.^p.n).^2;
        nm = sqrt(gF1.^2+gF9.^2)+1e-15;
        n1r=gF1./nm; n9r=gF9./nm;
        F1=(1/p.epsilon)*(p.vs*(1+p.mu)*Kn./(Kn+x9.^p.n)-p.vm*x1./(p.Km+x1));
        F9_idx=9;
        F9=p.k3_p*x(:,10)-p.k4_p*x9;
        M_simp(i)=trapz(t,w.*(F1.*n1r+F9.*n9r));

        % Lyapunov surrogate: L=integral w*<gradV, F> dt
        % V=0.5*||x-x*||^2, gradV=x-x*, F=full RHS
        gV=x-repmat(x_star,size(x,1),1);
        F_all=zeros(size(x));
        for k=1:length(t)
            F_all(k,:)=circ_ode(t(k),x(k,:)',p)';
        end
        M_lyap(i)=trapz(t, w.*sum(gV.*F_all,2));
    catch
        M_simp(i)=NaN; M_lyap(i)=NaN;
    end
    if mod(i,5)==0, fprintf('Progress %d/%d\n',i,length(mu_range)); end
end

M_use=M_lyap;
idx=find(diff(sign(M_use)),1);
if ~isempty(idx)
    mu_c=interp1(M_use(idx:idx+1),mu_range(idx:idx+1),0);
    fprintf('Critical mu_c = %.4f\n',mu_c);
    fprintf('Paper: 0.523, MATCONT: 0.556\n');
else
    mu_c=NaN;
    fprintf('No zero crossing. L range: [%.4f, %.4f]\n',...
        min(M_lyap(~isnan(M_lyap))),max(M_lyap(~isnan(M_lyap))));
end

results.mu_range=mu_range;results.M_simp=M_simp;results.M_lyap=M_lyap;
results.mu_critical=mu_c;results.lambda=lambda;
results.dimension=10;results.model_name='Circadian Rhythm';

figure; plot(mu_range,M_lyap,'r','LineWidth',2); hold on;
yline(0,'k--'); grid on;
xlabel('\mu'); ylabel('L(\mu) Lyapunov surrogate');
title('Model 5: Circadian 10D');
if ~isnan(mu_c), xline(mu_c,'b--','LineWidth',1.5); end
save('circadian_10D_results.mat','results');
end

function dxdt=circ_ode(~,x,p)
dxdt=zeros(10,1);
Kn=p.KI^p.n;
dxdt(1)=(1/p.epsilon)*(p.vs*(1+p.mu)*Kn/(Kn+x(9)^p.n)-p.vm*x(1)/(p.Km+x(1)));
dxdt(2)=p.ks*x(1)-p.vd*x(2)/(p.Kd+x(2));
dxdt(3)=(1/p.delta)*(p.k1*x(2)-p.k2*x(3)-p.v1*x(3)/(p.K1+x(3)));
dxdt(4)=p.k3*x(5)-p.k4*x(4);
dxdt(5)=(1/p.epsilon)*(p.k2*x(3)-p.k3*x(5));
dxdt(6)=p.vs_p*Kn/(Kn+x(4)^p.n)-p.vm_p*x(6)/(p.Km_p+x(6));
dxdt(7)=p.ks_p*x(6)-p.vd_p*x(7)/(p.Kd_p+x(7));
dxdt(8)=(1/p.delta)*(p.k1_p*x(7)-p.k2_p*x(8)-p.v1_p*x(8)/(p.K1_p+x(8)));
dxdt(9)=p.k3_p*x(10)-p.k4_p*x(9);
dxdt(10)=(1/p.epsilon)*(p.k2_p*x(8)-p.k3_p*x(10));
end

function lambda=compute_lambda_circ(x0,p)
h=1e-5; n=10; J=zeros(n);
for i=1:n
    xp=x0;xm=x0;xp(i)=xp(i)+h;xm(i)=xm(i)-h;
    fp=circ_ode(0,xp,p);fm=circ_ode(0,xm,p);
    J(:,i)=(fp-fm)/(2*h);
end
e=eig(J);e=real(e);e=e(e<0);
if isempty(e), lambda=50; else, lambda=-max(e); end
end
