function dydt = enzymatic_rhs(~, y)
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

%% enzymatic_rhs.m  — ODE for enzymatic feedback CRN
% Authors:  Sri Venkata Durga Sudarsan Madhyannapu, Pradheep Kumar S.
% Journal: CNSNS (Elsevier), Manuscript ID: CNSNS-D-26-00848

x  = y(1); y1 = y(2); z = y(3);
dx  = 1.5*x*(1-x) - 0.4*x*z;
dy1 = 0.2*x - 0.8*y1;
dz  = 0.5*y1 - 0.6*z;
dydt = [dx; dy1; dz];
end
