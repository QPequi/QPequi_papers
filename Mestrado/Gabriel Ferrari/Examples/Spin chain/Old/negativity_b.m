%% Cadeia com campo externo.

%% Variando B

Neg = zeros(1000,1);
BB = zeros(1000,1);

for B = 0.001:.001:1

BB(round(B.*1000)) = B;
    
T = 0.01;
J = 1;

n = B./T;
m=J./(2*T);
Z = 1 + 2.*cosh(2.*n) + 4.*cosh(n).*cosh(m) + 4.*cosh(sqrt(2).*m);

% Rho parcialmente transposto.

M = sinh(J./(2.*T));
N = cosh(J./(2.*T));
P = sinh(J./(sqrt(2).*T));
Q = cosh(J./(sqrt(2).*T));


rhot = (1./Z).* [exp(-2.*n) 0 0 0 exp(-n).*M 0 0 0 Q-(1./2); 0 exp(-n).*N 0 0 0 sqrt(2).*P 0 0 0; 0 0 Q+(1./2) 0 0 0 0 0 0; 0 0 0 exp(-n).*N 0 0 0 sqrt(2).*P 0; exp(-n).*M 0 0 0 2.*Q 0 0 0 exp(n).*M; 0 sqrt(2).*P 0 0 0 exp(n).*N 0 0 0; 0 0 0 0 0 0 Q+(1./2) 0 0; 0 0 0 sqrt(2).*P 0 0 0 exp(n).*N 0; Q-(1./2) 0 0 0 exp(n).*M 0 0 0 exp(2.*n)];

mod_rhot= trace(sqrt(rhot'*rhot));

% Negatividade

Neg(round(B.*1000)) = (mod_rhot-1)./2;

end

plot(BB,Neg)



%% Variando T 


% Neg = zeros(1000,1);
% TT = zeros(1000,1);
% 
% for T = 0.001:.001:1
% 
% TT(round(T.*1000)) = T;
%     
% B = 0;
% J = 1;
% 
% n = B./T;
% m=J./(2*T);
% Z = 1 + 2.*cosh(2.*n) + 4.*cosh(n).*cosh(m)+4.*cosh(sqrt(2).*m);
% 
% % Rho parcialmente transposto.
% 
% M = sinh(J./(2.*T));
% N = cosh(J./(2.*T));
% P = sinh(J./(sqrt(2).*T));
% Q = cosh(J./(sqrt(2).*T));
% 
% 
% rhot = (1./Z).* [exp(-2.*n) 0 0 0 exp(-n).*M 0 0 0 Q-(1./2); 0 exp(-n).*N 0 0 0 sqrt(2).*P 0 0 0; 0 0 Q+(1./2) 0 0 0 0 0 0; 0 0 0 exp(-n).*N 0 0 0 sqrt(2).*P 0; exp(-n).*M 0 0 0 2.*Q 0 0 0 exp(n).*M; 0 sqrt(2).*P 0 0 0 exp(n).*N 0 0 0; 0 0 0 0 0 0 Q+(1./2) 0 0; 0 0 0 sqrt(2).*P 0 0 0 exp(n).*N 0; Q-(1./2) 0 0 0 exp(n).*M 0 0 0 exp(2.*n)];
% 
% mod_rhot= trace(sqrt(rhot'*rhot));
% 
% % Negatividade
% 
% Neg(round(T.*1000)) = (mod_rhot-1)./2;
% 
% end
% 
% plot(TT,Neg)