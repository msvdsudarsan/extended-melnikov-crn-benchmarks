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
    
    % Initial condition on critical manifold
    x0 = [0.35, 0.28, 0.42, 0.31, 0.38, 0.29, 0.33];
    
    % Compute dominant eigenvalue for contraction rate
    lambda = compute_lambda_glycolysis(x0, params);
    fprintf('Dominant contraction rate lambda = %.2f\n', lambda);
    
    % Compute Melnikov functional for each mu
    for i = 1:length(mu_range)
        mu = mu_range(i);
        params.mu = mu;
        
        % Solve reduced system
        [t, x] = ode15s(@(t,x) glycolysis_reduced(t, x, params), tspan, x0, options);
        
        % Compute simplified Melnikov integral
        M_simp(i) = compute_melnikov_integral(t, x, params, lambda);
        
        if mod(i, 10) == 0
            fprintf('Progress: %d/%d\n', i, length(mu_range));
        end
    end
    
    % Find zero crossing
    zero_idx = find(diff(sign(M_simp)), 1);
    if ~isempty(zero_idx)
        mu_critical = interp1(M_simp(zero_idx:zero_idx+1), ...
                              mu_range(zero_idx:zero_idx+1), 0);
        fprintf('\nCritical parameter mu_c = %.4f\n', mu_critical);
    else
        mu_critical = NaN;
        fprintf('\nNo zero crossing found\n');
    end
    
    % Store results
    results.mu_range = mu_range;
    results.M_simp = M_simp;
    results.mu_critical = mu_critical;
    results.lambda = lambda;
    results.params = params;
    
    % Plot results
    figure('Position', [100 100 800 500]);
    plot(mu_range, M_simp, 'b-', 'LineWidth', 2);
    hold on;
    plot([mu_range(1) mu_range(end)], [0 0], 'k--', 'LineWidth', 1);
    if ~isnan(mu_critical)
        plot(mu_critical, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
        text(mu_critical, 0, sprintf(' \\mu_c = %.3f', mu_critical), ...
             'FontSize', 12, 'VerticalAlignment', 'bottom');
    end
    xlabel('\mu (perturbation parameter)', 'FontSize', 12);
    ylabel('M_{simp}(\mu)', 'FontSize', 12);
    title('Glycolysis 7D Model: Simplified Melnikov Function', 'FontSize', 14);
    grid on;
    
    % Save results
    save('glycolysis_7D_results.mat', 'results');
    saveas(gcf, 'glycolysis_7D_melnikov.png');
    fprintf('\nResults saved to glycolysis_7D_results.mat\n');
end

function dxdt = glycolysis_reduced(t, x, params)
    % 7-variable glycolysis model ODEs
    % x = [x1, x2, x3, x4, x5, x6, x7]
    
    v0 = params.v0;
    k1 = params.k1;
    k2 = params.k2;
    k3 = params.k3;
    k4 = params.k4;
    k5 = params.k5;
    k6 = params.k6;
    k7 = params.k7;
    k8 = params.k8;
    eps = params.epsilon;
    delta = params.delta;
    
    dxdt = zeros(7, 1);
    
    % Fast variable x1
    dxdt(1) = (1/eps) * (v0 - (k1 * x(1) * x(2)^2) / (1 + x(2)^2));
    
    % Slow variables
    dxdt(2) = (k1 * x(1) * x(2)^2) / (1 + x(2)^2) - k2 * x(2) - k3 * x(2) * x(3);
    dxdt(3) = k3 * x(2) * x(3) - k4 * x(3);
    
    % Medium-fast variable x4
    dxdt(4) = (1/delta) * (k4 * x(3) - k5 * x(4) * x(5));
    
    % Slow variables
    dxdt(5) = k5 * x(4) * x(5) - k6 * x(5);
    dxdt(6) = k6 * x(5) - k7 * x(6) * x(7);
    
    % Medium-fast variable x7
    dxdt(7) = (1/delta) * (k7 * x(6) * x(7) - k8 * x(7));
end

function lambda = compute_lambda_glycolysis(x0, params)
    % Compute dominant stable eigenvalue of fast subsystem Jacobian
    
    % Numerical Jacobian of fast dynamics
    eps = params.epsilon;
    h = 1e-6;
    
    % Fast subsystem includes x1, x4, x7 (fast variables)
    % Linearize around x0
    J = zeros(7, 7);
    
    for i = 1:7
        x_plus = x0;
        x_plus(i) = x_plus(i) + h;
        f_plus = glycolysis_reduced(0, x_plus, params);
        
        x_minus = x0;
        x_minus(i) = x_minus(i) - h;
        f_minus = glycolysis_reduced(0, x_minus, params);
        
        J(:, i) = (f_plus - f_minus) / (2*h);
    end
    
    % Compute eigenvalues
    eigs_J = eig(J);
    
    % Select dominant stable eigenvalue (most negative real part)
    stable_eigs = eigs_J(real(eigs_J) < 0);
    if isempty(stable_eigs)
        lambda = 5.0; % Default fallback
        warning('No stable eigenvalues found, using default lambda = 5.0');
    else
        [~, idx] = max(abs(real(stable_eigs)));
        lambda = -real(stable_eigs(idx));
    end
    
    fprintf('Eigenvalue analysis: lambda = %.4f\n', lambda);
end

function M = compute_melnikov_integral(t, x, params, lambda)
    % Compute simplified Melnikov integral with exponential weighting
    
    % Perturbation vector p = [x2, 0, 0, 0, 0, 0, 0]'
    p_vals = x(:, 2);
    
    % Compute normal vector (simplified for high-dimensional case)
    % Use gradient of slow manifold projection
    dx = diff(x(:, 1:3), 1, 1);
    dt_diff = diff(t);
    
    % Normalize
    n_vals = zeros(length(t)-1, 1);
    for i = 1:length(t)-1
        tangent = dx(i, :) / dt_diff(i);
        n_vals(i) = norm(tangent);
    end
    n_vals = [n_vals; n_vals(end)]; % Extend to match length
    
    % Exponential weight
    w = exp(-lambda * (t - t(1)));
    
    % Inner product <p, n>
    integrand = w .* p_vals .* n_vals;
    
    % Trapezoidal integration
    M = trapz(t, integrand);
end

% Run the analysis
fprintf('====================================\n');
fprintf('Glycolysis 7D Model - Melnikov Analysis\n');
fprintf('Extreme stiffness: epsilon = 0.0008\n');
fprintf('====================================\n\n');

results = glycolysis_7D_melnikov();

fprintf('\n====================================\n');
fprintf('Analysis complete!\n');
fprintf('====================================\n');
