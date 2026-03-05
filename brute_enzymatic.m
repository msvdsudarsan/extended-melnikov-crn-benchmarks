function R = brute_enzymatic()
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

%% brute_enzymatic.m
% High-accuracy reference simulation for enzymatic feedback error metric
%
% Paper: "An Efficient Melnikov-Lyapunov Computational Diagnostic for
%         Transition Detection in Slow-Fast Chemical Reaction Networks"
% Authors:  Sri Venkata Durga Sudarsan Madhyannapu, Pradheep Kumar S.
% Journal:  Communications in Nonlinear Science and Numerical Simulation
%          (Elsevier), ISSN: 1007-5704
% Manuscript ID: CNSNS-D-26-00848
% Status: Under Review, 2026

opts = odeset('RelTol',1e-9,'AbsTol',1e-10);
[~,X] = ode15s(@enzymatic_rhs, [0 30], [1.0 0.2 0.3], opts);
R     = max(abs(X(:,1)));
end
