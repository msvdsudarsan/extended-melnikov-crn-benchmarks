% Model 3: Extended Glycolysis Model (7 Variables)
% Extreme stiffness regime: epsilon = 0.0008
% Reference: Goldbeter, A. (1996) Biochemical Oscillations and Cellular Rhythms

function results = glycolysis_7D_melnikov()

    % Parameters
    params.v0 = 0.36;
    params.k1 = 1.0;
    params.k2 = 6.0;
    params.k3 = 0.8;
    params.k4 = 1.2;
    params.k5 = 0.9;
    params.k6 = 1.5;
    params.k7 = 0.7;
    params.k8 = 1.1;
    params.epsilon = 0.0008;  % Extreme stiffness
    params.delta = 0.005;

    % Perturbation parameter range
    mu_range = linspace(0.1, 0.5, 50);

    % Preallocate
    M_simp = zeros(size(mu_range));

    % Integration parameters
    tspan = [0, 40];
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);

    % Initial condition
    x0 = [0.35, 0.28, 0.42, 0.31, 0.38, 0.29, 0.33];

    % Compute contraction rate
    lambda = compute_lambda_glycolysis(x0, params);
    fprintf('Dominant contraction rate lambda = %.2f\n', lambda);

    % Loop over mu
    for i = 1:length(mu_range)
        params.mu = mu_range(i);

        [t, x] = ode15s(@(t,x) glycolysis_reduced(t, x, params), ...
                        tspan, x0, options);

        M_simp(i) = compute_melnikov_integral(t, x, lambda);

        if mod(i,10)==0
            fprintf('Progress: %d/%d\n', i, length(mu_range));
        end
    end

    % Zero crossing
    zero_idx = find(diff(sign(M_simp)), 1);
    if ~isempty(zero_idx)
        mu_critical = interp1(M_simp(zero_idx:zero_idx+1), ...
                              mu_range(zero_idx:zero_idx+1), 0);
    else
        mu_critical = NaN;
    end

    % Store results
    results.mu_range = mu_range;
    results.M_simp = M_simp;
    results.mu_critical = mu_critical;
    results.lambda = lambda;

    % Plot
    figure;
    plot(mu_range, M_simp,'b','LineWidth',2); hold on;
    yline(0,'k--');
    xlabel('\mu'); ylabel('M_{simp}');
    title('Glycolysis 7D Melnikov Function');
    grid on;

    save('glycolysis_7D_results.mat','results');
end

% ---------- ODE SYSTEM ----------
function dxdt = glycolysis_reduced(~, x, p)

dxdt = zeros(7,1);

dxdt(1) = (1/p.epsilon)*(p.v0 - (p.k1*x(1)*x(2)^2)/(1+x(2)^2));
dxdt(2) = (p.k1*x(1)*x(2)^2)/(1+x(2)^2) - p.k2*x(2) - p.k3*x(2)*x(3);
dxdt(3) = p.k3*x(2)*x(3) - p.k4*x(3);
dxdt(4) = (1/p.delta)*(p.k4*x(3) - p.k5*x(4)*x(5));
dxdt(5) = p.k5*x(4)*x(5) - p.k6*x(5);
dxdt(6) = p.k6*x(5) - p.k7*x(6)*x(7);
dxdt(7) = (1/p.delta)*(p.k7*x(6)*x(7) - p.k8*x(7));
end

% ---------- LAMBDA ----------
function lambda = compute_lambda_glycolysis(x0,p)

h=1e-6; J=zeros(7);
for i=1:7
 xp=x0; xm=x0;
 xp(i)=xp(i)+h; xm(i)=xm(i)-h;
 fp=glycolysis_reduced(0,xp,p);
 fm=glycolysis_reduced(0,xm,p);
 J(:,i)=(fp-fm)/(2*h);
end

e=eig(J);
e=e(real(e)<0);
if isempty(e)
 lambda=5;
else
 lambda=-max(real(e));
end
end

% ---------- MELNIKOV ----------
function M = compute_melnikov_integral(t,x,lambda)

p_vals=x(:,2);
dx=diff(x(:,1:3));
dt=diff(t);

n=zeros(length(t)-1,1);
for i=1:length(n)
 n(i)=norm(dx(i,:)/dt(i));
end
n=[n;n(end)];

w=exp(-lambda*(t-t(1)));
M=trapz(t, w.*p_vals.*n);
end
