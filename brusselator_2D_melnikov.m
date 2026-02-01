% Model 7: Brusselator (2D, Benchmark)
% Extreme stiffness regime: epsilon = 0.001
% Reference: Prigogine & Lefever (1968)

function results = brusselator_2D_melnikov()

    % Parameters
    params.A = 1.0;
    params.B = 3.0;
    params.epsilon = 0.001;

    mu_range = linspace(2.5, 3.5, 50);
    M_simp = zeros(size(mu_range));

    tspan = [0 30];
    opts  = odeset('RelTol',1e-9,'AbsTol',1e-11);
    x0    = [params.A, params.B/params.A];

    lambda = compute_lambda_brusselator(params);
    fprintf('Dominant contraction rate lambda = %.3f\n', lambda);

    for i = 1:length(mu_range)
        params.B = mu_range(i);

        [t,x] = ode15s(@(t,x) brusselator_reduced(t,x,params), tspan, x0, opts);

        M_simp(i) = compute_melnikov_integral_brusselator(t,x,lambda);

        if mod(i,10)==0
            fprintf('Progress %d/%d\n', i, length(mu_range));
        end
    end

    idx = find(diff(sign(M_simp)),1);
    if isempty(idx)
        mu_c = NaN;
    else
        mu_c = interp1(M_simp(idx:idx+1),mu_range(idx:idx+1),0);
    end

    results.mu_range = mu_range;
    results.M_simp   = M_simp;
    results.mu_critical = mu_c;
    results.lambda   = lambda;
    results.dimension = 2;
    results.model_name = 'Brusselator';

    figure;
    plot(mu_range,M_simp,'b','LineWidth',2); hold on;
    yline(0,'k--');
    xlabel('B'); ylabel('M_{simp}');
    title('Brusselator 2D Melnikov Function');
    grid on;

    save('brusselator_2D_results.mat','results');
end

% ================= ODE =================
function dxdt = brusselator_reduced(~,x,p)

dxdt = zeros(2,1);

dxdt(1) = (1/p.epsilon)*(p.A - (p.B+1)*x(1) + x(1)^2*x(2));
dxdt(2) = p.B*x(1) - x(1)^2*x(2);
end

% ================= LAMBDA =================
function lambda = compute_lambda_brusselator(p)

A=p.A; B=p.B; eps=p.epsilon;

x_eq=A; y_eq=B/A;

J=[-(B+1)+2*x_eq*y_eq, x_eq^2;
   B-2*x_eq*y_eq,      -x_eq^2];

J(1,:) = J(1,:)/eps;

e=eig(J);
e=e(real(e)<0);
if isempty(e)
 lambda=3.5;
else
 lambda=-max(real(e));
end
end

% ================= MELNIKOV =================
function M = compute_melnikov_integral_brusselator(t,x,lambda)

p_vals=x(:,1);
dx=diff(x); dt=diff(t);

n=zeros(length(t)-1,1);
for i=1:length(n)
 v=dx(i,:)/dt(i);
 n(i)=norm(v);
end
n=[n;n(end)];

w=exp(-lambda*(t-t(1)));
M=trapz(t, w.*p_vals.*n);
end
