% Model 5: Mammalian Circadian Rhythm Network (10 Variables)
% Extreme stiffness regime: epsilon = 0.0002
% Reference: Leloup & Goldbeter (2003) PNAS

function results = circadian_10D_melnikov()
    % Parameters - Circadian clock model
    params.vs = 1.6;
    params.vm = 0.8;
    params.Km = 0.5;
    params.ks = 0.6;
    params.vd = 0.95;
    params.Kd = 0.3;
    params.k1 = 0.4;
    params.k2 = 0.2;
    params.KI = 1.0;
    params.n = 4;  % Hill coefficient
    params.v1 = 0.7;
    params.K1 = 0.5;
    params.k3 = 0.4;
    params.k4 = 0.2;
    
    % Prime parameters (second gene)
    params.vs_p = 1.5;
    params.vm_p = 0.75;
    params.Km_p = 0.5;
    params.ks_p = 0.65;
    params.vd_p = 1.0;
    params.Kd_p = 0.35;
    params.k1_p = 0.45;
    params.k2_p = 0.22;
    params.v1_p = 0.75;
    params.K1_p = 0.55;
    params.k3_p = 0.42;
    params.k4_p = 0.21;
    
    params.epsilon = 0.0002;  % Extreme stiffness
    params.delta = 0.002;
    
    % Perturbation parameter range
    mu_range = linspace(0.3, 0.7, 40);
    
    % Preallocate
    M_simp = zeros(size(mu_range));
    M_lyap = zeros(size(mu_range)); % Lyapunov surrogate for comparison
    
    % Integration parameters
    tspan = [0, 60];
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-11, 'MaxStep', 0.1);
    
    % Initial condition on critical manifold
    x0 = [0.5, 0.4, 0.45, 0.35, 0.42, 0.48, 0.38, 0.44, 0.36, 0.40];
    
    % Compute dominant eigenvalue
    lambda = compute_lambda_circadian(x0, params);
    fprintf('Dominant contraction rate lambda = %.2f\n', lambda);
    
    % Equilibrium for Lyapunov surrogate
    x_eq = x0; % Use initial as reference
    
    % Compute Melnikov functionals
    fprintf('\nComputing Melnikov functionals (10D system - may take time)...\n');
    
    for i = 1:length(mu_range)
        mu = mu_range(i);
        params.mu = mu;
        
        try
            % Solve reduced system
            [t, x] = ode15s(@(t,x) circadian_reduced(t, x, params), tspan, x0, options);
            
            % Standard Melnikov integral
            M_simp(i) = compute_melnikov_integral_circadian(t, x, params, lambda);
            
            % Lyapunov surrogate for extreme stiffness
            M_lyap(i) = compute_lyapunov_surrogate(t, x, params, lambda, x_eq);
            
        catch ME
            warning('Integration failed at mu = %.3f: %s', mu, ME.message);
            M_simp(i) = NaN;
            M_lyap(i) = NaN;
        end
        
        if mod(i, 5) == 0
            fprintf('Progress: %d/%d (mu = %.3f)\n', i, length(mu_range), mu);
        end
    end
    
    % Find zero crossings
    % Use Lyapunov surrogate for extreme stiffness regime
    M_use = M_lyap; % More stable for epsilon = 0.0002
    
    zero_idx = find(diff(sign(M_use)), 1);
    if ~isempty(zero_idx) && ~any(isnan(M_use(zero_idx:zero_idx+1)))
        mu_critical = interp1(M_use(zero_idx:zero_idx+1), ...
                              mu_range(zero_idx:zero_idx+1), 0);
        fprintf('\n*** Critical parameter mu_c = %.4f (Lyapunov surrogate) ***\n', mu_critical);
    else
        mu_critical = NaN;
        fprintf('\n*** No zero crossing found ***\n');
    end
    
    % Store results
    results.mu_range = mu_range;
    results.M_simp = M_simp;
    results.M_lyap = M_lyap;
    results.mu_critical = mu_critical;
    results.lambda = lambda;
    results.params = params;
    results.dimension = 10;
    results.model_name = 'Circadian Rhythm';
    
    % Plot results
    figure('Position', [100 100 1000 700]);
    
    subplot(2,2,1);
    plot(mu_range, M_simp, 'b-', 'LineWidth', 2);
    hold on;
    plot([mu_range(1) mu_range(end)], [0 0], 'k--', 'LineWidth', 1);
    xlabel('\mu', 'FontSize', 11);
    ylabel('M_{simp}', 'FontSize', 11);
    title('Standard Melnikov Functional', 'FontSize', 12);
    grid on;
    
    subplot(2,2,2);
    plot(mu_range, M_lyap, 'r-', 'LineWidth', 2);
    hold on;
    plot([mu_range(1) mu_range(end)], [0 0], 'k--', 'LineWidth', 1);
    if ~isnan(mu_critical)
        plot(mu_critical, 0, 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
        text(mu_critical, 0, sprintf(' \\mu_c = %.3f', mu_critical), ...
             'FontSize', 12, 'FontWeight', 'bold');
    end
    xlabel('\mu', 'FontSize', 11);
    ylabel('\mathcal{L}', 'FontSize', 11);
    title('Melnikov-Lyapunov Surrogate', 'FontSize', 12);
    grid on;
    
    subplot(2,2,[3,4]);
    % Time series at critical parameter
    if ~isnan(mu_critical)
        params.mu = mu_critical;
        try
            [t_crit, x_crit] = ode15s(@(t,x) circadian_reduced(t, x, params), ...
                                       [0 100], x0, options);
            plot(t_crit, x_crit(:,1), 'b-', 'LineWidth', 1.5);
            hold on;
            plot(t_crit, x_crit(:,6), 'r-', 'LineWidth', 1.5);
            xlabel('Time', 'FontSize', 11);
            ylabel('Concentration', 'FontSize', 11);
            title(sprintf('Time Series at \\mu_c = %.3f', mu_critical), 'FontSize', 12);
            legend('mRNA_1', 'mRNA_2', 'Location', 'best');
            grid on;
        catch
            text(0.5, 0.5, 'Integration failed at critical parameter', ...
                 'Units', 'normalized', 'HorizontalAlignment', 'center');
        end
    end
    
    sgtitle('Circadian Rhythm 10D Model: Extreme Stiffness Regime', ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    % Save results
    save('circadian_10D_results.mat', 'results');
    saveas(gcf, 'circadian_10D_melnikov.png');
    fprintf('\nResults saved to circadian_10D_results.mat\n');
end

function dxdt = circadian_reduced(t, x, params)
    % 10-variable mammalian circadian clock model
    % Variables: x1=mRNA_Per, x2=PER_cyt, x3=PER_nuc_P, x4=PER_nuc,
    %            x5=PER_complex, x6=mRNA_Cry, x7=CRY_cyt, x8=CRY_nuc_P,
    %            x9=CRY_nuc, x10=CRY_complex
    
    vs = params.vs; vm = params.vm; Km = params.Km;
    ks = params.ks; vd = params.vd; Kd = params.Kd;
    k1 = params.k1; k2 = params.k2; KI = params.KI;
    n = params.n; v1 = params.v1; K1 = params.K1;
    k3 = params.k3; k4 = params.k4;
    
    vs_p = params.vs_p; vm_p = params.vm_p; Km_p = params.Km_p;
    ks_p = params.ks_p; vd_p = params.vd_p; Kd_p = params.Kd_p;
    k1_p = params.k1_p; k2_p = params.k2_p;
    v1_p = params.v1_p; K1_p = params.K1_p;
    k3_p = params.k3_p; k4_p = params.k4_p;
    
    eps = params.epsilon;
    delta = params.delta;
    
    dxdt = zeros(10, 1);
    
    % PER pathway
    % x1: mRNA_Per (fast)
    dxdt(1) = (1/eps) * (vs * KI^n / (KI^n + x(9)^n) - vm * x(1) / (Km + x(1)));
    
    % x2: PER_cyt (slow)
    dxdt(2) = ks * x(1) - vd * x(2) / (Kd + x(2));
    
    % x3: PER_nuc_P (medium-fast)
    dxdt(3) = (1/delta) * (k1 * x(2) - k2 * x(3) - v1 * x(3) / (K1 + x(3)));
    
    % x4: PER_nuc (slow)
    dxdt(4) = k3 * x(5) - k4 * x(4);
    
    % x5: PER_complex (fast)
    dxdt(5) = (1/eps) * (k2 * x(3) - k3 * x(5));
    
    % CRY pathway
    % x6: mRNA_Cry (slow)
    dxdt(6) = vs_p * KI^n / (KI^n + x(4)^n) - vm_p * x(6) / (Km_p + x(6));
    
    % x7: CRY_cyt (slow)
    dxdt(7) = ks_p * x(6) - vd_p * x(7) / (Kd_p + x(7));
    
    % x8: CRY_nuc_P (medium-fast)
    dxdt(8) = (1/delta) * (k1_p * x(7) - k2_p * x(8) - v1_p * x(8) / (K1_p + x(8)));
    
    % x9: CRY_nuc (slow)
    dxdt(9) = k3_p * x(10) - k4_p * x(9);
    
    % x10: CRY_complex (fast)
    dxdt(10) = (1/eps) * (k2_p * x(8) - k3_p * x(10));
end

function lambda = compute_lambda_circadian(x0, params)
    % Compute dominant eigenvalue - simplified for 10D
    
    % For extreme stiffness, use analytical estimate
    % Fast timescale ~ 1/epsilon
    lambda = 1.0 / params.epsilon; % Initial estimate
    
    % Refine via numerical Jacobian (subsampled)
    h = 1e-6;
    n = length(x0);
    
    % Sample fast variables only for efficiency
    fast_idx = [1, 5, 10]; % mRNA_Per, PER_complex, CRY_complex
    
    J_fast = zeros(length(fast_idx), length(fast_idx));
    
    for i = 1:length(fast_idx)
        idx = fast_idx(i);
        x_plus = x0;
        x_plus(idx) = x_plus(idx) + h;
        f_plus = circadian_reduced(0, x_plus, params);
        
        x_minus = x0;
        x_minus(idx) = x_minus(idx) - h;
        f_minus = circadian_reduced(0, x_minus, params);
        
        J_fast(:, i) = (f_plus(fast_idx) - f_minus(fast_idx)) / (2*h);
    end
    
    eigs_fast = eig(J_fast);
    stable_eigs = eigs_fast(real(eigs_fast) < 0);
    
    if ~isempty(stable_eigs)
        lambda = -min(real(stable_eigs)); % Use most stable
    end
    
    fprintf('Eigenvalue (fast subsystem): lambda = %.2f\n', lambda);
end

function M = compute_melnikov_integral_circadian(t, x, params, lambda)
    % Standard Melnikov integral
    
    % Perturbation on mRNA variables: p = [x1, 0, ..., 0]'
    p_vals = x(:, 1);
    
    % Normal estimate via fast variable derivatives
    dx_fast = [diff(x(:,1)); 0];
    n_vals = abs(dx_fast);
    n_vals = n_vals / max(n_vals + 1e-10);
    
    % Exponential weight
    w = exp(-lambda * (t - t(1)));
    
    % Integrand
    integrand = w .* p_vals .* (1 + n_vals);
    
    % Integration
    M = trapz(t, integrand);
end

function L = compute_lyapunov_surrogate(t, x, params, lambda, x_eq)
    % Melnikov-Lyapunov surrogate for extreme stiffness
    
    % Perturbation
    p_vals = x(:, 1);
    
    % Lyapunov gradient: ∇V = x - x_eq
    grad_V = x - repmat(x_eq, size(x, 1), 1);
    grad_V_norm = sqrt(sum(grad_V.^2, 2));
    
    % Exponential weight
    w = exp(-lambda * (t - t(1)));
    
    % Inner product <∇V, p>
    integrand = w .* p_vals .* grad_V_norm;
    
    % Integration
    L = trapz(t, integrand);
end

% Execute analysis
fprintf('================================================\n');
fprintf('CIRCADIAN RHYTHM 10D MODEL\n');
fprintf('Melnikov-Lyapunov Analysis\n');
fprintf('Extreme stiffness: epsilon = 0.0002\n');
fprintf('================================================\n\n');

tic;
results = circadian_10D_melnikov();
elapsed = toc;

fprintf('\n================================================\n');
fprintf('Analysis complete in %.2f seconds\n', elapsed);
fprintf('Dimension: %dD\n', results.dimension);
fprintf('Model: %s\n', results.model_name);
fprintf('Lyapunov surrogate used for extreme stiffness\n');
fprintf('================================================\n');
