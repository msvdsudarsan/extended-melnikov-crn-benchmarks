% Model 4: MAPK Cascade (8 Variables)
% Extreme stiffness regime: epsilon = 0.0005
% Reference: Kholodenko, B.N. (2000) Eur. J. Biochem.

function results = mapk_cascade_8D_melnikov()

    % Parameters
    params.k1 = 1.2;
    params.k2 = 0.8;
    params.k3 = 1.5;
    params.k4 = 0.9;
    params.k5 = 1.1;
    params.k6 = 0.7;
    params.k7 = 1.3;
    params.k8 = 0.85;
    params.k9 = 1.0;
    params.S  = 1.0;
    params.epsilon = 0.0005;
    params.delta   = 0.003;

    % Parameter sweep
    mu_range = linspace(0.2, 0.6, 50);
    M_simp   = zeros(size(mu_range));

    % Integration setup
    tspan = [0 50];
    opts  = odeset('RelTol',1e-9,'AbsTol',1e-11);
    x0    = [0.45 0.38 0.42 0.35 0.40 0.33 0.37 0.31];

    % Contraction rate
    lambda = compute_lambda_mapk(x0, params);
    fprintf('Dominant contraction rate lambda = %.3f\n', lambda);

    % Loop over mu
    for i = 1:length(mu_range)
        params.mu = mu_range(i);

        [t,x] = ode15s(@(t,x) mapk_reduced(t,x,params), tspan, x0, opts);

        M_simp(i) = compute_melnikov_integral_mapk(t,x,lambda);

        if mod(i,10)==0
            fprintf('Progress %d/%d\n',i,length(mu_range));
        end
    end

    % Zero crossing
    idx = find(diff(sign(M_simp)),1);
    if isempty(idx)
        mu_c = NaN;
    else
        mu_c = interp1(M_simp(idx:idx+1),mu_range(idx:idx+1),0);
    end

    % Results struct
    results.mu_range    = mu_range;
    results.M_simp      = M_simp;
    results.mu_critical = mu_c;
    results.lambda      = lambda;
    results.dimension   = 8;
    results.model_name  = 'MAPK Cascade';

    % Plot
    figure;
    plot(mu_range,M_simp,'b','LineWidth',2); hold on;
    yline(0,'k--');
    xlabel('\mu'); ylabel('M_{simp}');
    title('MAPK Cascade 8D Melnikov Function');
    grid on;

    save('mapk_8D_results.mat','results');
end

% ================= ODE =================
function dxdt = mapk_reduced(~,x,p)

dxdt = zeros(8,1);

dxdt(1) = (1/p.epsilon)*(p.k1*p.S - p.k2*x(1)*x(2));
dxdt(2) = p.k2*x(1)*x(2) - p.k3*x(2);

dxdt(3) = (1/p.epsilon)*(p.k3*x(2) - p.k4*x(3)*x(4));
dxdt(4) = p.k4*x(3)*x(4) - p.k5*x(4);

dxdt(5) = (1/p.epsilon)*(p.k5*x(4) - p.k6*x(5)*x(6));
dxdt(6) = p.k6*x(5)*x(6) - p.k7*x(6);

dxdt(7) = (1/p.delta)*(p.k7*x(6) - p.k8*x(7)*x(8));
dxdt(8) = p.k8*x(7)*x(8) - p.k9*x(8);
end

% ================= LAMBDA =================
function lambda = compute_lambda_mapk(x0,p)

h = 1e-7; n = length(x0); J = zeros(n);
for i=1:n
 xp=x0; xm=x0;
 xp(i)=xp(i)+h; xm(i)=xm(i)-h;
 fp=mapk_reduced(0,xp,p);
 fm=mapk_reduced(0,xm,p);
 J(:,i)=(fp-fm)/(2*h);
end

e = eig(J);
e = e(real(e)<0);
if isempty(e)
 lambda = 10;
else
 lambda = -max(real(e));
end
end

% ================= MELNIKOV =================
function M = compute_melnikov_integral_mapk(t,x,lambda)

p_vals = x(:,1);
dx = diff(x(:,[1 3 5]));
dt = diff(t);

n=zeros(length(t)-1,1);
for i=1:length(n)
 n(i)=norm(dx(i,:)/dt(i));
end
n=[n;n(end)];

w=exp(-lambda*(t-t(1)));
M=trapz(t, w.*p_vals.*n);
end
