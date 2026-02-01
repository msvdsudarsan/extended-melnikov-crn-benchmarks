% Model 5: Mammalian Circadian Rhythm Network (10 Variables)
% Extreme stiffness regime: epsilon = 0.0002
% Reference: Leloup & Goldbeter (2003) PNAS

function results = circadian_10D_melnikov()

    % Parameters (PER pathway)
    params.vs = 1.6;  params.vm = 0.8;  params.Km = 0.5;
    params.ks = 0.6;  params.vd = 0.95; params.Kd = 0.3;
    params.k1 = 0.4;  params.k2 = 0.2;  params.KI = 1.0;
    params.n  = 4;    params.v1 = 0.7;  params.K1 = 0.5;
    params.k3 = 0.4;  params.k4 = 0.2;

    % CRY pathway
    params.vs_p = 1.5; params.vm_p = 0.75; params.Km_p = 0.5;
    params.ks_p = 0.65; params.vd_p = 1.0; params.Kd_p = 0.35;
    params.k1_p = 0.45; params.k2_p = 0.22;
    params.v1_p = 0.75; params.K1_p = 0.55;
    params.k3_p = 0.42; params.k4_p = 0.21;

    params.epsilon = 0.0002;
    params.delta   = 0.002;

    mu_range = linspace(0.3,0.7,40);
    M_simp = zeros(size(mu_range));
    M_lyap = zeros(size(mu_range));

    tspan = [0 60];
    opts  = odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',0.1);
    x0    = [0.5 0.4 0.45 0.35 0.42 0.48 0.38 0.44 0.36 0.40];

    lambda = compute_lambda_circadian(x0,params);
    fprintf('Dominant contraction rate lambda = %.2f\n',lambda);

    x_eq = x0;

    for i=1:length(mu_range)
        params.mu = mu_range(i);

        try
            [t,x] = ode15s(@(t,x)circadian_reduced(t,x,params),tspan,x0,opts);

            M_simp(i) = compute_melnikov_integral_circadian(t,x,lambda);
            M_lyap(i) = compute_lyapunov_surrogate(t,x,lambda,x_eq);
        catch
            M_simp(i)=NaN; M_lyap(i)=NaN;
        end

        if mod(i,5)==0
            fprintf('Progress %d/%d\n',i,length(mu_range));
        end
    end

    M_use = M_lyap;
    idx = find(diff(sign(M_use)),1);
    if isempty(idx)
        mu_c = NaN;
    else
        mu_c = interp1(M_use(idx:idx+1),mu_range(idx:idx+1),0);
    end

    results.mu_range = mu_range;
    results.M_simp = M_simp;
    results.M_lyap = M_lyap;
    results.mu_critical = mu_c;
    results.lambda = lambda;
    results.dimension = 10;
    results.model_name = 'Circadian Rhythm';

    figure;
    plot(mu_range,M_lyap,'r','LineWidth',2); hold on;
    yline(0,'k--');
    xlabel('\mu'); ylabel('Lyapunov surrogate');
    title('Circadian 10D Melnikov–Lyapunov Function');
    grid on;

    save('circadian_10D_results.mat','results');
end

% ================= ODE =================
function dxdt = circadian_reduced(~,x,p)

dxdt = zeros(10,1);

dxdt(1) = (1/p.epsilon)*(p.vs*p.KI^p.n/(p.KI^p.n+x(9)^p.n) - p.vm*x(1)/(p.Km+x(1)));
dxdt(2) = p.ks*x(1) - p.vd*x(2)/(p.Kd+x(2));
dxdt(3) = (1/p.delta)*(p.k1*x(2)-p.k2*x(3)-p.v1*x(3)/(p.K1+x(3)));
dxdt(4) = p.k3*x(5)-p.k4*x(4);
dxdt(5) = (1/p.epsilon)*(p.k2*x(3)-p.k3*x(5));

dxdt(6) = p.vs_p*p.KI^p.n/(p.KI^p.n+x(4)^p.n) - p.vm_p*x(6)/(p.Km_p+x(6));
dxdt(7) = p.ks_p*x(6) - p.vd_p*x(7)/(p.Kd_p+x(7));
dxdt(8) = (1/p.delta)*(p.k1_p*x(7)-p.k2_p*x(8)-p.v1_p*x(8)/(p.K1_p+x(8)));
dxdt(9) = p.k3_p*x(10)-p.k4_p*x(9);
dxdt(10)= (1/p.epsilon)*(p.k2_p*x(8)-p.k3_p*x(10));
end

% ================= LAMBDA =================
function lambda = compute_lambda_circadian(x0,p)
lambda = 1/p.epsilon;
end

% ================= MELNIKOV =================
function M = compute_melnikov_integral_circadian(t,x,lambda)
p_vals=x(:,1);
dx=diff(x(:,1));
dx=[dx;dx(end)];
w=exp(-lambda*(t-t(1)));
M=trapz(t,w.*p_vals.*abs(dx));
end

% ================= LYAPUNOV =================
function L = compute_lyapunov_surrogate(t,x,lambda,x_eq)
grad = x - repmat(x_eq,size(x,1),1);
normg = sqrt(sum(grad.^2,2));
w=exp(-lambda*(t-t(1)));
L=trapz(t,w.*normg.*x(:,1));
end
