%% 
clear all
clc

%% Geradores G_ij e F_ij

% Geradores do Grupo SU(3)


Ger = cell(8,1);
for j = 1:8
    Gen{j} = zeros(3,3);
end
    


Gen{1} = (1./2).*[0 1 0;1 0 0; 0 0 0];
Gen{2} = (1./2).*[0 -1i 0; 1i 0 0; 0 0 0];
Gen{3} = (1./2).*[1 0 0; 0 -1 0; 0 0 0];
Gen{4} = (1./2).*[0 0 1; 0 0 0; 1 0 0];
Gen{5} = (1./2).*[0 0 -1i; 0 0 0; 1i 0 0];
Gen{6} = (1./2).*[0 0 0; 0 0 1; 0 1 0];
Gen{7} = (1./2).*[0 0 0; 0 0 -1i; 0 1i 0];
Gen{8} = (1./(2.*sqrt(3))).*[1 0 0; 0 1 0; 0 0 -2];
II = eye(3,3);


% g(ijk)

for i = 1:8
    for j = 1:8
        for k = 1:8
            g(i,j,k) = (1./4).*trace((Gen{i}*Gen{j}+Gen{j}*Gen{i})*Gen{k});
        end
    end
end

% for i = 1:8
%     for j = 1:8
%         for k = 1:8
%             f(i,j,k) = (1./4.*1i).*trace((Gen{i}*Gen{j}-Gen{j}*Gen{i})*Gen{k});
%         end
%     end
% end
% 



G = cell(8,8);


for i = 1:8
    for j = 1:8
        G{i,j} = [g(i,j,1) g(i,j,2) g(i,j,3) g(i,j,4) g(i,j,5) g(i,j,6) g(i,j,7) g(i,j,8)];
    end
end


%% Definição de B

B = 1



%% Matrizz W:
W = zeros(8,8);

passo1 = 0.01;
Ni1 = 0.01;
Nf1 = 2.01;
D1 = (Nf1 - Ni1)./passo1
passo2 = .1;
Ni2 = -10.0;
Nf2 = 10.0;
D2 = (Nf2 - Ni2)./passo2

TT = zeros(201,1);
JJ = zeros(201,1);
W = zeros(8,8);
LQU = zeros(201,201);

%  for T =Ni1:passo1:Nf1
%         for J = Ni2:passo2:Nf2
%            
%             TT(round(T.*100)) = T;
%             JJ(round(J.*10 + 101)) = J;
%             
%         end
%  end
% 


raiz_rho = zeros(9,9);
% Calculo das correlacoes:
    for T =Ni1:passo1:Nf1
        for J = Ni2:passo2:Nf2
           
            TT(round(T.*100)) = T;
            JJ(round(J.*10 + 101)) = J;
            
%            Calculo da raiz de rho:
                rhoo_b;
                raiz_rho = sqrt(rho);
                
                L = [trace(rho*kron(Gen{1},II)); trace(rho*kron(Gen{2},II)); trace(rho*kron(Gen{3},II)); trace(rho*kron(Gen{4},II)); trace(rho*kron(Gen{5},II)); trace(rho*kron(Gen{6},II)); trace(rho*kron(Gen{7},II)); trace(rho*kron(Gen{8},II))];
                
            for y=1:8
                for z=1:8
                    W(y,z) = trace(raiz_rho*kron(Gen{y},eye(3))*raiz_rho*kron(Gen{z},eye(3))) - G{y,z}*L;
                end
            end
            
             U = (2./3)-max(real(eig(W)));
             LQU(round(T.*100),round(J.*10 + 101)) = U;
           
        end    
    end        


mesh(JJ,TT,LQU)




clear B D1 D2 Gen Ger II J L Nf1 Nf2 Ni1 Ni2 T U d f g i j k passo1 passo2 raiz_rho rho u y z
