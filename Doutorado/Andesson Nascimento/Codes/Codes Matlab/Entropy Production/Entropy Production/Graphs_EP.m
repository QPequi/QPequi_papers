clear all
cla

load('vect.mat');
load('EP_j300_h02.mat');
load('EP_j300_h08.mat');
load('EP_j300_h04.mat');
load('EP_j300_h06.mat');
load('LB_EP_j300_h02.mat');
load('LB_EP_j300_h08.mat');
load('UB_EP_j300_h02.mat');
load('UB_EP_j300_h08.mat');


% Defina o Interpreter globalmente
set(groot, 'defaultLegendInterpreter', 'latex');

%% Plot dos gráficos

plot(vect,EP_j300_h02,'Linewidth',2.5,'LineStyle','--','DisplayName','$h=0,2$','Color', 'r')
hold on
plot(vect,EP_j300_h04,'Linewidth',2.5,'LineStyle','-.','DisplayName','$h=0,4$','Color', 'g')
hold on
plot(vect,EP_j300_h06,'Linewidth',2.5,'LineStyle','-','DisplayName','$h=0,6$','Color', 'k')
hold on
plot(vect,EP_j300_h08,'Linewidth',2.5,'LineStyle','-','DisplayName','$h=0,8$','Color', 'b')
hold off
ylim([0,40])
xlim([-1,30])
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$s(\frac{2}{\pi}\mathcal{L})$', 'Interpreter', 'latex')
ax = gca;
legend('Location', 'best');
ax.FontSize = 18; % Altere para o tamanho desejado
ax.Box = 'on'








