% Código utilizado para plotar os gráficos do Angulo de Bures
% calculado para o modelo LMG (código 'Bures.m')
clear all
clc

load('vect.mat');
load('beta.mat');
load('BA_j100_beta_h02.mat');
load('BA_j100_beta_h08.mat');
load('BA_j300_beta_h02.mat');
BA_j300_beta_h02 = LL;
clear LL
load('BA_j300_beta_h08.mat');
load('BA_j300_beta001_h02.mat');
load('BA_j300_beta001_h08.mat');

% Defina o Interpreter globalmente
set(groot, 'defaultLegendInterpreter', 'latex');

%% Plot do Angulo de Bures para quench h = 0,2 -- j=100

plot(vect, BA_j100_beta_h02(1,:), 'LineWidth', 2,'LineStyle',':','Color','b','DisplayName', '$\beta=0,1$');
hold on
plot(vect, BA_j100_beta_h02(2,:), 'LineWidth', 2,'LineStyle','-.','Color','r','DisplayName', '$\beta=0,5$');
hold on
plot(vect, BA_j100_beta_h02(3,:), 'LineWidth', 2,'LineStyle',':','Color','k','DisplayName', '$\beta=1$');
hold on
plot(vect, BA_j100_beta_h02(4,:), 'LineWidth', 2,'LineStyle',':','Color','b','DisplayName', '$\beta=2$');
hold on
plot(vect, BA_j100_beta_h02(5,:), 'LineWidth', 2,'LineStyle','--','Color','k','DisplayName', '$\beta=4$');
hold on
plot(vect, BA_j100_beta_h02(6,:), 'LineWidth', 2,'LineStyle',':','Color','m','DisplayName', '$\beta=6$');
hold off
ylim([0.7,1.6])
xlim([-1,20])
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\mathcal{L}(\rho_\tau,\rho_{\tau}^{eq})$', 'Interpreter', 'latex')
ax = gca;
legend('Location', 'best');
ax.FontSize = 18; % Altere para o tamanho desejado
ax.Box = 'on'

%% Plot do Angulo de Bures para quench h = 0,8 -- j=100

plot(vect, BA_j100_beta_h08(1,:), 'LineWidth', 2,'LineStyle','-','Color','k','DisplayName', '$\beta=0,1$');
hold on
plot(vect, BA_j100_beta_h08(2,:), 'LineWidth', 2,'LineStyle','-.','Color','r','DisplayName', '$\beta=0,5$');
hold on
plot(vect, BA_j100_beta_h08(3,:), 'LineWidth', 2,'LineStyle',':','Color','g','DisplayName', '$\beta=1$');
hold on
plot(vect, BA_j100_beta_h08(4,:), 'LineWidth', 2,'LineStyle',':','Color','b','DisplayName', '$\beta=2$');
hold on
plot(vect, BA_j100_beta_h08(5,:), 'LineWidth', 2,'LineStyle','--','Color','k','DisplayName', '$\beta=4$');
hold on
plot(vect, BA_j100_beta_h08(6,:), 'LineWidth', 2,'LineStyle',':','Color','m','DisplayName', '$\beta=6$');
hold off
ylim([0.7,1.6])
xlim([-1,20])
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\mathcal{L}(\rho_\tau,\rho_{\tau}^{eq})$', 'Interpreter', 'latex')
ax = gca;
legend('Location', 'best');
ax.FontSize = 18; % Altere para o tamanho desejado
ax.Box = 'on'

%% Plot do Angulo de Bures para quench h = 0,2 -- j=300

plot(vect, BA_j300_beta001_h02_, 'LineWidth', 2,'LineStyle','-.','Color','r','DisplayName', '$\beta=0,01$');
hold on
plot(vect, BA_j300_beta_h02(1,:), 'LineWidth', 2,'LineStyle','--','Color','b','DisplayName', '$\beta=0,1$');
hold on
% plot(vect, BA_j300_beta_h02(2,:), 'LineWidth', 2,'LineStyle','-.','Color','r','DisplayName', '$\beta=0,5$');
% hold on
plot(vect, BA_j300_beta_h02(3,:), 'LineWidth', 2,'LineStyle','-','Color','k','DisplayName', '$\beta=1$');
% hold on
% plot(vect, BA_j300_beta_h02(4,:), 'LineWidth', 2,'LineStyle',':','Color','b','DisplayName', '$\beta=2$');
%hold on
%plot(vect, BA_j300_beta_h02(5,:), 'LineWidth', 2,'LineStyle','--','Color','k','DisplayName', '$\beta=4$');
%hold on
%plot(vect, BA_j300_beta_h02(6,:), 'LineWidth', 2,'LineStyle',':','Color','m','DisplayName', '$\beta=6$');
hold off
ylim([0.9,1.6])
xlim([-1,30])
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\mathcal{L}(\rho_\tau,\rho_{\tau}^{eq})$', 'Interpreter', 'latex')
ax = gca;
legend('Location', 'best');
ax.FontSize = 18; % Altere para o tamanho desejado
ax.Box = 'on'

%% Plot do Angulo de Bures para quench h = 0,8 -- j=300

plot(vect, BA_j300_beta001_h08, 'LineWidth', 2,'LineStyle','-.','Color','r','DisplayName', '$\beta=0,01$');
hold on
plot(vect, BA_j300_beta_h08(1,:), 'LineWidth', 2,'LineStyle','--','Color','b','DisplayName', '$\beta=0,1$');
% hold on
% plot(vect, BA_j300_beta_h08(2,:), 'LineWidth', 2,'LineStyle','-.','Color','r','DisplayName', '$\beta=0,5$');
hold on
plot(vect, BA_j300_beta_h08(3,:), 'LineWidth', 2,'LineStyle','-','Color','k','DisplayName', '$\beta=1$');
% hold on
% plot(vect, BA_j300_beta_h08(4,:), 'LineWidth', 2,'LineStyle',':','Color','b','DisplayName', '$\beta=2$');
%hold on
%plot(vect, BA_j300_beta_h08(5,:), 'LineWidth', 2,'LineStyle','--','Color','k','DisplayName', '$\beta=4$');
%hold on
%plot(vect, BA_j00_beta_h08(6,:), 'LineWidth', 2,'LineStyle',':','Color','m','DisplayName', '$\beta=6$');
hold off
ylim([0.9,1.6])
xlim([-1,30])
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\mathcal{L}(\rho_\tau,\rho_{\tau}^{eq})$', 'Interpreter', 'latex')
ax = gca;
legend('Location', 'best');
ax.FontSize = 18; % Altere para o tamanho desejado
ax.Box = 'on'




