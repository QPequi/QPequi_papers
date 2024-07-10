clear all
clc

dimh = 50;
vech = linspace(0,1,dimh);

load('TAh50_j500_t1000_10.mat');
load('TAh50_j200_t1000_10.mat');
load('TAh50_j100_t1000_10.mat');

%% Plot

% Posição da linha vertical
x_line = 0.5;
% Intervalo do eixo y (ajuste de acordo com a necessidade)
y_range = [0.5,3.3];
x_range = [0,1];

% Defina o Interpreter globalmente
set(groot, 'defaultLegendInterpreter', 'latex');
% Plotar a linha vertical
figure;
hold on;
plot([x_line, x_line], y_range, 'g-.', 'LineWidth', 1.5,'DisplayName', '$h_c^d$');

plot(vech,TAh50_j100_t1000_10,'LineWidth',2,'DisplayName', '$j = 50$','Color','b','LineStyle',':')
hold on;
plot(vech,TAh50_j200_t1000_10,'LineWidth',2,'DisplayName', '$j = 100$','Color','r','LineStyle','--')
hold on;
plot(vech,TA_j500_t1000_10,'LineWidth',2,'DisplayName', '$j = 500$', 'Color','k','LineStyle','-')
legend('Location', 'best');
hold off;
xlabel('$h$', 'Interpreter', 'latex','FontSize',16);
ylabel('$\overline{\langle \Sigma \rangle}$', 'Interpreter', 'latex','FontSize',16)
ax = gca;
ax.FontSize = 16; % Altere para o tamanho desejado
% Defina os limites do eixo x e y
xlim(x_range);
ylim(y_range);
ax.Box = 'on'





