clear all
clc
disp('Running');

dimx = 100;
vecx = linspace(0,1,dimx);

for x=1:length(vecx)
% Calculo do limite inferior e superior
         s_lb(x) = 2*vecx(x)^2;
         s_ub(x) = - log(1-vecx(x));
               
        % Calculo de s(x) exato
      
        e = 1*10^(-4);
        vecr = linspace(vecx(x)+e,1-e,dimx);
        for r=1:length(vecr)
           s_exact(x,r) = (vecr(r)-vecx(x))*log((vecr(r)-vecx(x))/vecr(r)) + ...
           + [1-vecr(r)+vecx(x)]*log((1-vecr(r)+vecx(x))/(1-vecr(r)));
        end
    end
    for x=1:length(vecx)
        s_exact_min(x)=min(s_exact(x,:));
    end   

% Defina o Interpreter globalmente
set(groot, 'defaultLegendInterpreter', 'latex');
    
plot(vecx,s_exact_min,'LineWidth',3,'LineStyle','-','Color','k');
hold on;
plot(vecx,s_lb,'LineWidth',3,'LineStyle',':','Color','b');
hold on;
plot(vecx,s_ub,'Linewidth',3,'LineStyle','--','Color','r');
hold off;
xlabel('$x$','Interpreter', 'latex');
ylabel('$s(x)$', 'Interpreter', 'latex')
ax = gca;
ax.FontSize = 20; % Altere para o tamanho desejado
%legend('Location', 'best');
ylim([-0.05,4.5])
xlim([-0.01,1.01])

disp('Ok');

