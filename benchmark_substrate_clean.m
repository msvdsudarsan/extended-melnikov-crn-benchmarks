function benchmark_substrate_clean()
% Benchmark for substrate inhibition oscillator — Extended Melnikov
disp("=== Substrate inhibition benchmark ===");

tic;
M = melnikov_substrate();
tM = toc;

tic;
R = brute_substrate();
tR = toc;

err = abs(M - R);

fprintf("Extended Melnikov runtime = %.6f sec\n", tM);
fprintf("Error metric              = %.6e\n", err);
end
