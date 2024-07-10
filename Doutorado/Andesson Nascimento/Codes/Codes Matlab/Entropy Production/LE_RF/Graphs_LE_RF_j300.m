% Código utilizado para plotar os gráficos do Eco de Loschmiidt e da função
% taxa calculados para o modelo LMG (código 'Loshmidt_Echo.m')
clear all
clc

load('vect.mat');
load('LE_j300_h02.mat');
load('RF_j300_h02.mat');
load('LE_j300_h03.mat');
load('RF_j300_h03.mat');
load('LE_j300_h08.mat');
load('RF_j300_h08.mat');

% Defina o Interpreter globalmente
set(groot, 'defaultLegendInterpreter', 'latex');

%% Plot do Eco de Loschmidt

plot(vect, LE_j300_h02, 'LineWidth', 2,'LineStyle','-.','Color','r','DisplayName', '$h=0,2$');
hold on
plot(vect, LE_j300_h08, 'LineWidth', 2,'LineStyle','-','Color','b','DisplayName', '$h=0,8$');
hold off
ylim([-0.05,1.1])
xlim([-1,30])
xlabel('$t$', 'Interpreter', 'latex','FontSize',14);
ylabel('$L_e$', 'Interpreter', 'latex')
ax = gca;
legend('Location', 'best');
ax.FontSize = 18; % Altere para o tamanho desejado
ax.Box = 'on'

%,'FontSize',28

%% Plot da Função Taxa

plot(vect, RF_j300_h02, 'LineWidth', 2,'LineStyle','-.','Color','r','DisplayName', '$h=0,2$');
hold on
plot(vect, RF_j300_h08, 'LineWidth', 2,'LineStyle','-','Color','b','DisplayName', '$h=0,8$');
%hold off
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\lambda(t)$', 'Interpreter', 'latex')
ylim([-0.005,0.13])
xlim([-1,30])
ax = gca;
legend('Location', 'best');
ax.FontSize = 18; % Altere para o tamanho desejado
ax.Box = 'on'










