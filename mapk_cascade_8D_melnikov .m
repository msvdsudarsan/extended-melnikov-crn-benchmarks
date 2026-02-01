% Model 4: MAPK Cascade (8 Variables)
% Extreme stiffness regime: epsilon = 0.0005
% Reference: Kholodenko, B.N. (2000) Eur. J. Biochem.

function results = mapk_cascade_8D_melnikov()
    % Parameters - MAPK cascade kinetic constants
    params.k1 = 1.2;
    params.k2 = 0.8;
    params.k3 = 1.5;
    params.k4 = 0.9;
    params.k5 = 1.1;
    params.k6 = 0.7;
    params.k7 = 1.3;
    params.k8 = 0.85;
    params.k9 = 1.0;
    params.S = 1.0;  % Stimulus strength
    params.epsilon = 0.0005;  % Extreme stiffness
    params.delta = 0.003;
    
    % Perturbation parameter range
    mu_range = linspace(0.2, 0.6, 50);
    
    % Preallocate
    M_simp = zeros(size(mu_range));
    
    % Integration parameters
    tspan = [0, 50];
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
    
    % Initial condition on critical manifold
    x0 = [0.45, 0.38, 0.42, 0.35, 0.40, 0.33, 0.37, 0.31];
    
    % Compute dominant eigenvalue for contraction rate
    lambda = compute_lambda_mapk(x0, params);
    fprintf('Dominant contraction rate lambda = %.2f\n', lambda);
    
    % Compute Melnikov functional for each mu
    fprintf('\nComputing Melnikov function across parameter range...\n');
    for i = 1:length(mu_range)
        mu = mu_range(i);
        params.mu = mu;
        
        % Solve reduced system
        [t, x] = ode15s(@(t,x) mapk_reduced(t, x, params), tspan, x0, options);
        
        % Compute simplified Melnikov integral
        M_simp(i) = compute_melnikov_integral_mapk(t, x, params, lambda);
        
        if mod(i, 10) == 0
            fprintf('Progress: %d/%d (mu = %.3f, M = %.4e)\n', ...
                    i, length(mu_range), mu, M_simp(i));
        end
    end
    
    % Find zero crossing
    zero_idx = find(diff(sign(M_simp)), 1);
    if ~isempty(zero_idx)
        mu_critical = interp1(M_simp(zero_idx:zero_idx+1), ...
                              mu_range(zero_idx:zero_idx+1), 0);
        fprintf('\n*** Critical parameter mu_c = %.4f ***\n', mu_critical);
    else
        mu_critical = NaN;
        fprintf('\n*** No zero crossing found ***\n');
    end
    
    % Store results
    results.mu_range = mu_range;
    results.M_simp = M_simp;
    results.mu_critical = mu_critical;
    results.lambda = lambda;
    results.params = params;
    results.dimension = 8;
    results.model_name = 'MAPK Cascade';
    
    % Plot results
    figure('Position', [100 100 900 600]);
    
    subplot(2,1,1);
    plot(mu_range, M_simp, 'b-', 'LineWidth', 2.5);
    hold on;
    plot([mu_range(1) mu_range(end)], [0 0], 'k--', 'LineWidth', 1.5);
    if ~isnan(mu_critical)
        plot(mu_critical, 0, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
        text(mu_critical, min(M_simp)*0.1, sprintf(' \\mu_c = %.3f', mu_critical), ...
             'FontSize', 13, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
    end
    xlabel('\mu (perturbation parameter)', 'FontSize', 12);
    ylabel('M_{simp}(\mu)', 'FontSize', 12);
    title('MAPK Cascade 8D: Simplified Melnikov Function', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    subplot(2,1,2);
    % Solve at critical parameter for phase portrait
    if ~isnan(mu_critical)
        params.mu = mu_critical;
        [t_crit, x_crit] = ode15s(@(t,x) mapk_reduced(t, x, params), [0 100], x0, options);
        plot(x_crit(:,1), x_crit(:,2), 'r-', 'LineWidth', 1.5);
        xlabel('x_1 (MAPKKK)', 'FontSize', 11);
        ylabel('x_2 (MAPKKK-P)', 'FontSize', 11);
        title(sprintf('Phase Portrait at Critical Parameter (\\mu = %.3f)', mu_critical), ...
              'FontSize', 12);
        grid on;
    end
    
    % Save results
    save('mapk_8D_results.mat', 'results');
    saveas(gcf, 'mapk_8D_melnikov.png');
    fprintf('\nResults saved to mapk_8D_results.mat\n');
    fprintf('Figure saved to mapk_8D_melnikov.png\n');
end

function dxdt = mapk_reduced(t, x, params)
    % 8-variable MAPK cascade ODEs
    % x1, x3, x5 = fast variables (phosphorylated forms at each tier)
    % x2, x4, x6, x8 = slow variables
    % x7 = medium-fast variable
    
    k1 = params.k1; k2 = params.k2; k3 = params.k3;
    k4 = params.k4; k5 = params.k5; k6 = params.k6;
    k7 = params.k7; k8 = params.k8; k9 = params.k9;
    S = params.S;
    eps = params.epsilon;
    delta = params.delta;
    
    dxdt = zeros(8, 1);
    
    % Tier 1: MAPKKK activation (fast)
    dxdt(1) = (1/eps) * (k1 * S - k2 * x(1) * x(2));
    dxdt(2) = k2 * x(1) * x(2) - k3 * x(2);
    
    % Tier 2: MAPKK activation (fast)
    dxdt(3) = (1/eps) * (k3 * x(2) - k4 * x(3) * x(4));
    dxdt(4) = k4 * x(3) * x(4) - k5 * x(4);
    
    % Tier 3: MAPK activation (fast)
    dxdt(5) = (1/eps) * (k5 * x(4) - k6 * x(5) * x(6));
    dxdt(6) = k6 * x(5) * x(6) - k7 * x(6);
    
    % Feedback regulation (medium-fast)
    dxdt(7) = (1/delta) * (k7 * x(6) - k8 * x(7) * x(8));
    dxdt(8) = k8 * x(7) * x(8) - k9 * x(8);
end

function lambda = compute_lambda_mapk(x0, params)
    % Compute dominant stable eigenvalue via numerical Jacobian
    
    h = 1e-7;
    n = length(x0);
    J = zeros(n, n);
    
    % Numerical differentiation
    for i = 1:n
        x_plus = x0;
        x_plus(i) = x_plus(i) + h;
        f_plus = mapk_reduced(0, x_plus, params);
        
        x_minus = x0;
        x_minus(i) = x_minus(i) - h;
        f_minus = mapk_reduced(0, x_minus, params);
        
        J(:, i) = (f_plus - f_minus) / (2*h);
    end
    
    % Eigenvalue analysis
    eigs_J = eig(J);
    
    % Select dominant stable eigenvalue
    stable_eigs = eigs_J(real(eigs_J) < -0.01); % Small threshold to avoid near-zero
    if isempty(stable_eigs)
        lambda = 10.0; % Default for stiff systems
        warning('No clearly stable eigenvalues, using default lambda = 10.0');
    else
        [~, idx] = max(abs(real(stable_eigs)));
        lambda = -real(stable_eigs(idx));
    end
    
    fprintf('Eigenvalue analysis complete: lambda = %.4f\n', lambda);
    fprintf('Number of stable eigenvalues: %d\n', length(stable_eigs));
end

function M = compute_melnikov_integral_mapk(t, x, params, lambda)
    % Compute simplified Melnikov integral
    
    % Perturbation: p = [x1, 0, 0, 0, 0, 0, 0, 0]'
    p_vals = x(:, 1);
    
    % Normal vector computation (using velocity magnitude as proxy)
    dx_fast = diff(x(:, [1,3,5]), 1, 1); % Fast variables
    dt_diff = diff(t);
    
    velocity_norm = zeros(length(t)-1, 1);
    for i = 1:length(t)-1
        velocity_norm(i) = norm(dx_fast(i, :)) / dt_diff(i);
    end
    velocity_norm = [velocity_norm; velocity_norm(end)];
    
    % Normalize to avoid numerical issues
    velocity_norm = velocity_norm / max(velocity_norm);
    
    % Exponential weight
    w = exp(-lambda * (t - t(1)));
    
    % Integrand
    integrand = w .* p_vals .* (1 + velocity_norm);
    
    % Trapezoidal integration
    M = trapz(t, integrand);
end

% Execute analysis
fprintf('========================================\n');
fprintf('MAPK CASCADE 8D MODEL\n');
fprintf('Melnikov Analysis for Bifurcation Detection\n');
fprintf('Extreme stiffness: epsilon = 0.0005\n');
fprintf('========================================\n\n');

tic;
results = mapk_cascade_8D_melnikov();
elapsed = toc;

fprintf('\n========================================\n');
fprintf('Analysis complete in %.2f seconds\n', elapsed);
fprintf('Dimension: %dD\n', results.dimension);
fprintf('Model: %s\n', results.model_name);
fprintf('========================================\n');
