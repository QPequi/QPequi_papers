clear all
clc

% Código para calcular Eco de Loschmidt e Função Taxa
dimt = 300;

vect = linspace(0,30,dimt);

% Definição das matrizes de Pauli
sigma_x = [0 1; 1 0];
sigma_y = [0 -1i; 1i 0];
sigma_z = [1 0; 0 -1];

g = 0.5; % Definição de g 
jj = 300; % Valor de j (número de spins)
dim = 2*jj + 1; % Dimensão da matriz

for m1 = -jj:jj
    for m2 = -jj:jj
        p1 = m1 + jj + 1;
        p2 = m2 + jj + 1;
        
        if m1 == m2 + 1
            matriz_x(p1, p2) = 0.5 * sqrt(jj*(jj+1)-m2*(m2+1));
        end
        if m1 == m2 - 1
            matriz_x(p1, p2) = 0.5 * sqrt(jj*(jj+1)-m2*(m2-1));
        end
        
        if m1 == m2 + 1
            matriz_y(p1, p2) = 0.5i * sqrt(jj*(jj+1)-m2*(m2+1));
        end
        if m1 == m2 - 1
            matriz_y(p1, p2) = -0.5i * sqrt(jj*(jj+1)-m2*(m2-1));
        end
        
        if m1 == m2
            matriz_z(p1, p2) = -m1;
        end
    end
clear m1
clear m2
clear p1
clear p2

end

% Hamiltoniano Inicial
h0 = 0; % definição do h inicial
H0 = - (1/jj)*((matriz_x*matriz_x) + g*(matriz_y*matriz_y)) - 2*h0*matriz_z;

% Autovalores e Autovetores
[V0, E0] = eig(H0);

% Hamiltoniano Final
h = 0.2;
H = -(1/jj)*((matriz_x*matriz_x) + g*(matriz_y*matriz_y)) - 2*h*matriz_z;
for t=1:length(vect)
    t
    U = expm(-1i*H*vect(t));
    L(t) = abs(V0(:,1)'*U*V0(:,1))^2; 
    RF(t) = -(1/(2*jj))*log(L(t));
end

%% Plot

plot(vect,L,'b','LineWidth',2.5);
%plot(vect,RF,'r','LineWidth',2.5);
xlabel('t','FontSize',18);
ylabel('$L_e$', 'Interpreter', 'latex','FontSize',16)
ax = gca;
ax.FontSize = 16;
%ylim([-0.05,1.1])
%xlim([-1,15])



disp('Cálculo concluído.');


