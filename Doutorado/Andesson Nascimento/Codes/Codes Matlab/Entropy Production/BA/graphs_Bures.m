% Código utilizado para plotar os gráficos do Angulo de Bures
% calculado para o modelo LMG (código 'Bures.m')
clear all
clc

load('vect.mat');
load('BA_j300_h02.mat');
load('BA_j300_h04.mat');
load('BA_j300_h06.mat');
load('BA_j300_h08.mat');

% Defina o Interpreter globalmente
set(groot, 'defaultLegendInterpreter', 'latex');

% Plot do Eco de Loschmidt

plot(vect, BA_j300_h02, 'LineWidth', 2,'LineStyle','-','Color','k','DisplayName', '$h=0,2$');
hold on
plot(vect, BA_j300_h04, 'LineWidth', 2,'LineStyle','-.','Color','r','DisplayName', '$h=0,4$');
hold on
plot(vect, BA_j300_h06, 'LineWidth', 2,'LineStyle',':','Color','g','DisplayName', '$h=0,6$');
hold on
plot(vect, BA_j300_h08, 'LineWidth', 2,'LineStyle','--','Color','b','DisplayName', '$h=0,8$');
hold off
ylim([0.9,1.6])
xlim([-1,30])
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\mathcal{L}(\rho_\tau,\rho_{\tau}^{eq})$', 'Interpreter', 'latex')
ax = gca;
legend('Location', 'best');
ax.FontSize = 18; % Altere para o tamanho desejado
ax.Box = 'on'