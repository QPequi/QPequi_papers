clear all
clc
disp('Running');

dimt = 300;
vect = linspace(0,30,dimt);

g = 0.5; % Definição de g 
jj = 300; % Valor de j (número de spins)
dim = 2*jj + 1; % Dimensão da matriz

% Inicialização das matrizes
matriz_x = zeros(dim);
matriz_y = zeros(dim);
matriz_z = zeros(dim);

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
end

% Hamiltoniano Inicial
h0 = 0;
H0 = - (1/jj)*((matriz_x*matriz_x) + g*(matriz_y*matriz_y)) - 2*h0*matriz_z;

% Autovalores e Autovetores
[V0, E0] = eig(H0);

% Estado térmico inicial
beta = [0.1,0.5,1,2,4,6];
for b =1:length(beta)
    tic
    beta(b)
    aux0 = expm(-beta(b)*H0);
    Z0 = trace(aux0);
    rho_th0 = aux0/Z0;
     % Substituir valores NaN por zero
    min_val = realmin;
    rho_th0(isnan(rho_th0)) = min_val;
    sq_rho0 = sqrtm(rho_th0);

%     % Construir a matriz rho_th0
%     rho_th0 = V0 * diag(rho_th0_diag) * V0';
    h = 0.8;
    H = -(1/jj)*((matriz_x*matriz_x) + g*(matriz_y*matriz_y)) - 2*h*matriz_z; 
    for t=1:length(vect)
        t
        U = expm(-1i*H*vect(t));
        psi_t = U*V0(:,1);
        rho_th = psi_t*psi_t';
        sq_rho_prod = sqrtm(sq_rho0*rho_th*sq_rho0);
        min_val2 = realmin;
        sq_rho_prod(isnan(sq_rho_prod)) = min_val2;
       
    
        % Raiz da Fidelidade
        sq_F = sqrtm(trace(sq_rho_prod)*trace(sq_rho_prod));
        LL(b,t) = acos(sq_F);
    end 
    toc
end
plot(vect,LL,'Linewidth',2)%,'DisplayName','exact')
set(gca, 'FontSize', 14);
%legend('Location', 'best');
xlabel('t','FontSize',16);
ylabel('$\mathcal{L}$', 'Interpreter', 'latex','FontSize',16)

disp('Ok');

clear sigma_x
clear sigma_y
clear sigma_z