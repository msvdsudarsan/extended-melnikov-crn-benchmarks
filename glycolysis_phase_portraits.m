clear; clc; close all;

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
