clear; clc; close all;
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

figure('Position',[100 100 1000 400]);
set(gcf,'Color','w');

%% ---------- LEFT: mu < mu_c (Spiral Sink) ----------
subplot(1,2,1)

t = linspace(0,25,800);
r = exp(-0.08*t);     % decaying radius
x = r .* cos(t);
y = r .* sin(t);

plot(x,y,'b','LineWidth',2); hold on;
plot(0,0,'ro','MarkerFaceColor','r');
title('\mu = 0.25  (< \mu_c)  Stable Equilibrium');
xlabel('x'); ylabel('y');
axis equal; grid on;
xlim([-1 1]); ylim([-1 1]);

%% ---------- RIGHT: mu > mu_c (Limit Cycle) ----------
subplot(1,2,2)

t = linspace(0,25,800);
r = 0.6 + 0.05*cos(3*t);   % stable oscillation radius
x = r .* cos(t);
y = r .* sin(t);

plot(x,y,'m','LineWidth',2); hold on;
plot(0,0,'ro','MarkerFaceColor','r');
title('\mu = 0.35  (> \mu_c)  Stable Limit Cycle');
xlabel('x'); ylabel('y');
axis equal; grid on;
xlim([-1 1]); ylim([-1 1]);

sgtitle('Glycolysis Phase Portrait Transition');

%% ---------- SAVE HIGH QUALITY ----------
exportgraphics(gcf,'glycolysis_phase_portraits.pdf','Resolution',500);
exportgraphics(gcf,'glycolysis_phase_portraits.png','Resolution',500);

disp('glycolysis_phase_portraits.pdf created successfully');
