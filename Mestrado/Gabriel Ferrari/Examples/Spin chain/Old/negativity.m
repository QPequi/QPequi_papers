%% Cadeia sem campo externo

% Variando J
% 
% passo = 0.0001;
% Nf = 10.0;
% Ni = -10.0;
% N = (Nf - Ni)/passo;
% Neg = zeros(N,1);
% JJ = zeros(N,1);
% 
% T = 100.0;
% k = 0;
% for J = Ni:passo:Nf
% 
% k = k + 1;
% JJ(k) = J;    
% 
% m=J./(2*T);
% Z = 3.0 + 4.*cosh(m) + 4.*cosh(sqrt(2).*m);
% 
% Rho parcialmente transposto.
% 
% M = sinh(m);
% N = cosh(m);
% P = sinh(sqrt(2).*m);
% Q = cosh(sqrt(2).*m);
% 
% 
% rhot = (1.0./Z).* [1 0 0 0 M 0 0 0 Q-0.5; 0 N 0 0 0 sqrt(2).*P 0 0 0; 0 0 Q+0.5 0 0 0 0 0 0; 0 0 0 N 0 0 0 sqrt(2).*P 0; M 0 0 0 2.*Q 0 0 0 M; 0 sqrt(2).*P 0 0 0 N 0 0 0; 0 0 0 0 0 0 Q+0.5 0 0; 0 0 0 sqrt(2).*P 0 0 0 N 0; Q-0.5 0 0 0 M 0 0 0 1];
% 
% mod_rhot= trace(sqrt((rhot')*(rhot)));
% 
% Negatividade
% 
% Neg(k) = (mod_rhot-1.0)./2.0;
% 
% end
% 
% plot(JJ,Neg)

%% Variando T

passo = 0.0001;
Nf = 2.0;
Ni = 0.0;
N = (Nf - Ni)/passo;
Neg = zeros(N,1);
TT = zeros(N,1);

k = 0;
J = 0.5;

for T = Ni:passo:Nf

k = k + 1;    
TT(k)= T;

m=J./(2*T);
Z = 3 + 4.*cosh(m) + 4.*cosh(sqrt(2).*m);

% Rho parcialmente transposto.

M = sinh(J./(2.*T));
N = cosh(J./(2.*T));
P = sinh(J./(sqrt(2).*T));
Q = cosh(J./(sqrt(2).*T));


rhot = (1./Z).* [1 0 0 0 M 0 0 0 Q-0.5; 0 N 0 0 0 sqrt(2).*P 0 0 0; 0 0 Q+0.5 0 0 0 0 0 0; 0 0 0 N 0 0 0 sqrt(2).*P 0; M 0 0 0 2.*Q 0 0 0 M; 0 sqrt(2).*P 0 0 0 N 0 0 0; 0 0 0 0 0 0 Q+0.5 0 0; 0 0 0 sqrt(2).*P 0 0 0 N 0; Q-0.5 0 0 0 M 0 0 0 1];

mod_rhot= trace(sqrt(rhot'*rhot));

% Negatividade

Neg(k) = (mod_rhot-1)./2;


end

plot(TT,Neg)

%%
