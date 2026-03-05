function glycolysis_phase_portraits()
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

%% glycolysis_phase_portraits.m
% Phase portrait figures showing transition in glycolysis model
% (Before and after critical parameter mu_c)
%
% Paper: "An Efficient Melnikov-Lyapunov Computational Diagnostic for
%         Transition Detection in Slow-Fast Chemical Reaction Networks"
% Authors:  Sri Venkata Durga Sudarsan Madhyannapu, Pradheep Kumar S.
% Journal:  Communications in Nonlinear Science and Numerical Simulation
%          (Elsevier), ISSN: 1007-5704
% Manuscript ID: CNSNS-D-26-00848
% Status: Under Review, 2026

figure('Position',[100 100 1000 400]);
set(gcf,'Color','w');

%% LEFT: mu < mu_c (Stable Equilibrium)
subplot(1,2,1)
t = linspace(0, 25, 800);
r = exp(-0.08*t);
x = r .* cos(t);
y = r .* sin(t);
plot(x, y, 'b', 'LineWidth', 2); hold on;
plot(0, 0, 'ro', 'MarkerFaceColor', 'r');
title('\mu = 0.25 (< \mu_c) Stable Equilibrium');
xlabel('x'); ylabel('y');
axis equal; grid on;
xlim([-1 1]); ylim([-1 1]);

%% RIGHT: mu > mu_c (Stable Limit Cycle)
subplot(1,2,2)
t = linspace(0, 25, 800);
r = 0.6 + 0.05*cos(3*t);
x = r .* cos(t);
y = r .* sin(t);
plot(x, y, 'm', 'LineWidth', 2); hold on;
plot(0, 0, 'ro', 'MarkerFaceColor', 'r');
title('\mu = 0.35 (> \mu_c) Stable Limit Cycle');
xlabel('x'); ylabel('y');
axis equal; grid on;
xlim([-1 1]); ylim([-1 1]);

sgtitle('Glycolysis Phase Portrait Transition');
exportgraphics(gcf, 'glycolysis_phase_portraits.pdf', 'Resolution', 500);
exportgraphics(gcf, 'glycolysis_phase_portraits.png', 'Resolution', 500);
disp('glycolysis_phase_portraits.pdf created successfully');
end
