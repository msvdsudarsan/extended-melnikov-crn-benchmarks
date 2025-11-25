function benchmark_enzymatic_clean()
% Benchmark for enzymatic feedback CRN
disp("=== Enzymatic feedback benchmark ===");

tic;
M = melnikov_enzymatic();
tM = toc;

tic;
R = brute_enzymatic();
tR = toc;

err = abs(M - R);

fprintf("Extended Melnikov runtime = %.6f sec\n", tM);
fprintf("Error metric              = %.6e\n", err);
end
