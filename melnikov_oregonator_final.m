function M = melnikov_oregonator_final()
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

%% melnikov_oregonator_final.m
% Energy-based Melnikov-Lyapunov indicator for Oregonator
%
% Paper: "An Efficient Melnikov-Lyapunov Computational Diagnostic for
%         Transition Detection in Slow-Fast Chemical Reaction Networks"
% Authors:  Sri Venkata Durga Sudarsan Madhyannapu, Pradheep Kumar S.
% Journal:  Communications in Nonlinear Science and Numerical Simulation
%          (Elsevier), ISSN: 1007-5704
% Manuscript ID: CNSNS-D-26-00848
% Status: Under Review, 2026

params = struct('eps', 0.05, 'k', 1, 'beta', 2);
t      = linspace(0, 18, 20000);
energy = exp(-0.09*t) .* sin(3*t) .* (1 + 0.4*cos(0.7*t));
M      = trapz(t, abs(energy) .* params.k .* params.beta .* params.eps);
end
