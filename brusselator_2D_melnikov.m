% Model 7: Brusselator (2D, Benchmark)
% Extreme stiffness regime: epsilon = 0.001
% Reference: Prigogine & Lefever (1968)

function results = brusselator_2D_melnikov()
    % Parameters
    params.A = 1.0;
    params.B = 3.0;
    params.epsilon = 0.001;  % Extreme stiffness
    
    % Perturbation parameter range
    mu_range = linspace(2.5, 3.5, 50);
    
    % Preallocate
    M_simp = zeros(size(mu_range));
    
    % Integration parameters
    tspan = [0, 30];
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
    
    % Initial condition on critical manifold
    x0 = [params.A, params.B/params.A];
    
    % Compute dominant eigenvalue for contraction rate
    lambda = compute_lambda_brusselator(x0, params);
    fprintf('Dominant contraction rate lambda = %.2f\n', lambda);
    
    % Compute Melnikov functional for each mu
    fprintf('\nComputing Melnikov function across parameter range...\n');
    for i = 1:length(mu_range)
        mu = mu_range(i);
        params.mu = mu;
        params.B = mu;  % B is the bifurcation parameter
        
        % Solve reduced system
        [t, x] = ode15s(@(t,x) brusselator_reduced(t, x, params), tspan, x0, options);
        
        % Compute simplified Melnikov integral
        M_simp(i) = compute_melnikov_integral_brusselator(t, x, params, lambda);
        
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
    results.dimension = 2;
    results.model_name = 'Brusselator';
    
    % Plot results
    figure('Position', [100 100 900 600]);
    
    subplot(2,1,1);
    plot(mu_range, M_simp, 'b-', 'LineWidth', 2.5);
    hold on;
    plot([mu_range(1) mu_range(end)], [0 0], 'k--', 'LineWidth', 1.5);
    if ~isnan(mu_critical)
        plot(mu_critical, 0, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
        text(mu_critical, min(M_simp)*0.1, sprintf(' B_c = %.3f', mu_critical), ...
             'FontSize', 13, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
    end
    xlabel('B (bifurcation parameter)', 'FontSize', 12);
    ylabel('M_{simp}(B)', 'FontSize', 12);
    title('Brusselator 2D: Simplified Melnikov Function', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    
    subplot(2,1,2);
    % Phase portrait at critical parameter
    if ~isnan(mu_critical)
        params.B = mu_critical;
        [t_crit, x_crit] = ode15s(@(t,x) brusselator_reduced(t, x, params), [0 50], x0, options);
        plot(x_crit(:,1), x_crit(:,2), 'r-', 'LineWidth', 1.5);
        xlabel('x (concentration)', 'FontSize', 11);
        ylabel('y (concentration)', 'FontSize', 11);
        title(sprintf('Phase Portrait at Critical Parameter (B = %.3f)', mu_critical), ...
              'FontSize', 12);
        grid on;
    end
    
    % Save results
    save('brusselator_2D_results.mat', 'results');
    saveas(gcf, 'brusselator_2D_melnikov.png');
    fprintf('\nResults saved to brusselator_2D_results.mat\n');
    fprintf('Figure saved to brusselator_2D_melnikov.png\n');
end

function dxdt = brusselator_reduced(t, x, params)
    % 2-variable Brusselator ODEs
    % x(1) = x, x(2) = y
    
    A = params.A;
    B = params.B;
    eps = params.epsilon;
    
    dxdt = zeros(2, 1);
    
    % Fast variable x (autocatalytic)
    dxdt(1) = (1/eps) * (A - (B + 1)*x(1) + x(1)^2 * x(2));
    
    % Slow variable y
    dxdt(2) = B*x(1) - x(1)^2 * x(2);
end

function lambda = compute_lambda_brusselator(x0, params)
    % Compute dominant stable eigenvalue
    
    A = params.A;
    B = params.B;
    eps = params.epsilon;
    
    % Jacobian of the Brusselator at equilibrium
    % x_eq = A, y_eq = B/A
    x_eq = A;
    y_eq = B/A;
    
    % Jacobian matrix
    J = [-(B+1) + 2*x_eq*y_eq, x_eq^2;
         B - 2*x_eq*y_eq, -x_eq^2];
    
    % Scale by 1/eps for fast variable
    J(1,:) = J(1,:) / eps;
    
    % Eigenvalues
    eigs_J = eig(J);
    
    % Select dominant stable eigenvalue
    stable_eigs = eigs_J(real(eigs_J) < -0.01);
    if isempty(stable_eigs)
        lambda = 3.5; % Default for Brusselator
        warning('No clearly stable eigenvalues, using default lambda = 3.5');
    else
        [~, idx] = max(abs(real(stable_eigs)));
        lambda = -real(stable_eigs(idx));
    end
    
    fprintf('Eigenvalue analysis complete: lambda = %.4f\n', lambda);
    fprintf('Eigenvalues: %.2f + %.2fi, %.2f + %.2fi\n', ...
            real(eigs_J(1)), imag(eigs_J(1)), real(eigs_J(2)), imag(eigs_J(2)));
end

function M = compute_melnikov_integral_brusselator(t, x, params, lambda)
    % Compute simplified Melnikov integral
    
    % Perturbation: p = [x, 0]'
    p_vals = x(:, 1);
    
    % Normal vector computation (perpendicular to flow)
    dx = diff(x, 1, 1);
    dt_diff = diff(t);
    
    % Velocity components
    v_x = zeros(length(t)-1, 1);
    v_y = zeros(length(t)-1, 1);
    for i = 1:length(t)-1
        v_x(i) = dx(i, 1) / dt_diff(i);
        v_y(i) = dx(i, 2) / dt_diff(i);
    end
    
    % Extend to match length
    v_x = [v_x; v_x(end)];
    v_y = [v_y; v_y(end)];
    
    % Normal vector (perpendicular): n = [-vy, vx] / norm
    n_x = -v_y;
    n_y = v_x;
    norm_n = sqrt(n_x.^2 + n_y.^2);
    
    % Avoid division by zero
    norm_n(norm_n < 1e-10) = 1;
    
    n_x = n_x ./ norm_n;
    
    % Exponential weight
    w = exp(-lambda * (t - t(1)));
    
    % Inner product <p, n> where p = [x, 0]'
    integrand = w .* p_vals .* n_x;
    
    % Trapezoidal integration
    M = trapz(t, integrand);
end

% Execute analysis
fprintf('========================================\n');
fprintf('BRUSSELATOR 2D MODEL\n');
fprintf('Melnikov Analysis - Benchmark Case\n');
fprintf('Extreme stiffness: epsilon = 0.001\n');
fprintf('========================================\n\n');

tic;
results = brusselator_2D_melnikov();
elapsed = toc;

fprintf('\n========================================\n');
fprintf('Analysis complete in %.2f seconds\n', elapsed);
fprintf('Dimension: %dD\n', results.dimension);
fprintf('Model: %s\n', results.model_name);
fprintf('========================================\n');
