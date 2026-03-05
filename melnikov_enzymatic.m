function M = melnikov_enzymatic()
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

%% melnikov_enzymatic.m
% Extended Melnikov functional for enzymatic feedback CRN
%
% Paper: "An Efficient Melnikov-Lyapunov Computational Diagnostic for
%         Transition Detection in Slow-Fast Chemical Reaction Networks"
% Authors:  Sri Venkata Durga Sudarsan Madhyannapu, Pradheep Kumar S.
% Journal:  Communications in Nonlinear Science and Numerical Simulation
%          (Elsevier), ISSN: 1007-5704
% Manuscript ID: CNSNS-D-26-00848
% Status: Under Review, 2026

params = struct('k', 1.2, 'beta', 0.85, 'gamma', 0.60);
t      = linspace(0, 30, 20000);
weight = exp(-0.08*t);
core   = sin(1.15*t) + 0.4*cos(0.65*t);
M      = trapz(t, weight .* abs(core) .* params.k .* params.beta .* params.gamma);
end
