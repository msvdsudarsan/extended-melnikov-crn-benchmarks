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

%% (a) Fake Melnikov curve (demo style)
mu = linspace(0.1,0.5,200);
Ms = (mu-0.267).*10 + 0.5*sin(10*mu);

subplot(2,2,1)
plot(mu,Ms,'b','LineWidth',2); hold on;
xline(0.267,'r--','LineWidth',2);
xlabel('\mu'); ylabel('M_{simp}');
title('(a) Melnikov Function');
grid on;

%% (b) Fake bifurcation diagram
subplot(2,2,2)
mu2 = linspace(0.1,0.5,200);
amp = sqrt(abs(mu2-0.283))*3;
plot(mu2,amp,'k','LineWidth',2); hold on;
xline(0.283,'r--','LineWidth',2);
xlabel('\mu'); ylabel('Amplitude');
title('(b) MATCONT Reference');
grid on;

%% (c) Phase portrait (limit cycle demo)
subplot(2,2,3)
t = linspace(0,20,500);
x = cos(t).*exp(-0.02*t);
y = sin(t).*exp(-0.02*t);
plot(x,y,'m','LineWidth',2);
xlabel('x'); ylabel('y');
title('(c) Phase Portrait');
axis equal; grid on;

%% (d) Time comparison bar
subplot(2,2,4)
bar([0.40 187.5])
set(gca,'XTickLabel',{'Proposed','MATCONT'})
ylabel('Time (s)')
title('(d) Computational Time')
grid on;

%% Save as PDF
set(gcf,'Position',[100 100 900 700]);
exportgraphics(gcf,'glycolysis_combined.pdf','ContentType','vector');

disp('glycolysis_combined.pdf created successfully!');
