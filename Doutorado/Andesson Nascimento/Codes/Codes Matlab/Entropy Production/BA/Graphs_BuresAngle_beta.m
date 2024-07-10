clear all
clc

load('vect.mat');
load('vecbeta.mat');
load('BA_j300_h02_beta.mat');
load('BA_j300_h08_beta.mat');

% Defina o Interpreter globalmente
set(groot, 'defaultLegendInterpreter', 'latex');

%% Plot do Angulo de Bures para quench h = 0,2

plot(vect, BA_j300_h02_beta(1,:), 'LineWidth', 2,'LineStyle','-','Color','k','DisplayName', '$\beta=0,1$');
hold on
plot(vect, BA_j300_h02_beta(5,:), 'LineWidth', 2,'LineStyle','-.','Color','r','DisplayName', '$\beta=0,5$');
hold on
plot(vect, BA_j300_h02_beta(10,:), 'LineWidth', 2,'LineStyle',':','Color','g','DisplayName', '$\beta=1$');
hold off
ylim([0.9,1.6])
xlim([-1,30])
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\mathcal{L}(\rho_\tau,\rho_{\tau}^{eq})$', 'Interpreter', 'latex')
ax = gca;
legend('Location', 'best');
ax.FontSize = 18; % Altere para o tamanho desejado
ax.Box = 'on'

%% Plot do Angulo de Bures para quench h = 0,8

plot(vect, BA_j300_h08_beta(1,:), 'LineWidth', 2,'LineStyle','-','Color','k','DisplayName', '$\beta=0,1$');
hold on
plot(vect, BA_j300_h08_beta(5,:), 'LineWidth', 2,'LineStyle','-.','Color','r','DisplayName', '$\beta=0,5$');
hold on
plot(vect, BA_j300_h08_beta(10,:), 'LineWidth', 2,'LineStyle',':','Color','g','DisplayName', '$\beta=1$');
hold off
ylim([0.9,1.6])
xlim([-1,30])
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\mathcal{L}(\rho_\tau,\rho_{\tau}^{eq})$', 'Interpreter', 'latex')
ax = gca;
legend('Location', 'best');
ax.FontSize = 18; % Altere para o tamanho desejado
ax.Box = 'on'
