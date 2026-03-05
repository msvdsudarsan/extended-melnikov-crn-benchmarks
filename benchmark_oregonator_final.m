function benchmark_oregonator_final()
%% Paper Title: "An Efficient Adjoint-Free Melnikov-Lyapunov Diagnostic for
%%               Transition Detection in Slow-Fast Chemical Reaction Networks"
%% Author 1:    Sri Venkata Durga Sudarsan Madhyannapu
%% Author 2:    Pradheep Kumar S.
%%
%% Affiliation 1: Freshmen Engineering Department, NRI Institute of Technology
%%                (Autonomous), Pothavarappadu, Agiripalli, Vijayawada,
%%                521212, Andhra Pradesh, India
%% Affiliation 2: Research Scholar, Jawaharlal Nehru Technological University
%%                Kakinada, Andhra Pradesh, India
%% Affiliation 3: School of Basic Sciences, SRM University AP, Neerukonda,
%%                Mangalagiri Mandal, Guntur, 522240, Andhra Pradesh, India
%%
%% Journal:       Communications in Nonlinear Science and Numerical Simulation
%%                Elsevier, ISSN: 1007-5704
%% Manuscript ID: CNSNS-D-26-00848
%% Status:        Under Review, 2026

%% benchmark_oregonator_final.m
% Benchmark for Oregonator — Extended Melnikov vs reference runtime
%
% Paper: "An Efficient Melnikov-Lyapunov Computational Diagnostic for
%         Transition Detection in Slow-Fast Chemical Reaction Networks"
% Authors:  Sri Venkata Durga Sudarsan Madhyannapu, Pradheep Kumar S.
% Journal:  Communications in Nonlinear Science and Numerical Simulation
%          (Elsevier), ISSN: 1007-5704
% Manuscript ID: CNSNS-D-26-00848
% Status: Under Review, 2026

disp('=== Balanced Oregonator benchmark ===');
tic;  M  = melnikov_oregonator_final();    tM = toc;
tic;  tR = brute_oregonator_final();
fprintf('Extended Melnikov runtime = %.6f sec\n', tM);
fprintf('Reference runtime         = %.6f sec\n', tR);
fprintf('Error metric              = Not applicable (energy-based)\n');
end
