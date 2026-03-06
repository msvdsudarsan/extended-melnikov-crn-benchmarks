function tR = brute_oregonator_final()
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
%% Submitted:      19 February 2026
%% Status:        Under Review, 2026
%% SSRN:          https://ssrn.com/abstract=6275667
%% SSRN ID:        6275667 (Distributed: 02/20/2026)
%%
% Reference runtime only (state-based error not meaningful)
opts = odeset("RelTol",1e-9,"AbsTol",1e-10);
tic;
ode15s(@oregonator_rhs,[0 2],[1 1 1],opts);
tR = toc;
end
