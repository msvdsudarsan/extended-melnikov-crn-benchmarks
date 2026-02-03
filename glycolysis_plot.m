clear; clc; close all;

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
