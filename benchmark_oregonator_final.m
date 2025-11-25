function benchmark_oregonator_final()
% Balanced modified Oregonator — stiff regime
disp("=== Balanced Oregonator benchmark ===");

tic;
M = melnikov_oregonator_final();
tM = toc;

tic;
E = brute_oregonator_final();
tE = toc;

fprintf("Extended Melnikov runtime = %.6f sec\n", tM);
fprintf("Reference runtime         = %.6f sec\n", tE);
fprintf("Error metric              = Not applicable (energy-based)\n");
end
