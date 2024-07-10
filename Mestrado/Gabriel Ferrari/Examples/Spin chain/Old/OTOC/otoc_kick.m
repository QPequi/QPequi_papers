% OTOC measure for the kicked top. The output of the program is the
% plot of the OTOC with respect to the kicks (time) for a random 
% initial condition.

% Global parameters 
p = pi/2; % Energy of the kicked top
N = 200; % Total number of kicks
tran = 1; % Transient to be disregarded

[sigmax,sigmay,sigmaz] = pauli(1/2); % Pauli matrices in z basis
j = 3/2; % Angular momemtum
dim = 2*j + 1; % Dimension of Hilbert space
  
[Jx,Jy,Jz] = pauli(j); % Angular momentum matrices in z basis

% Initial state vector
init = zeros(dim,1);
init(1) = 1;

% OTOC operators
epsilon = 0.0;
O1 = expm(-1i*p*Jz); % To be evolved
O2 = expm(-1i*(p+epsilon)*Jz); 

vec_kappa = 0 : 0.5 : 5;
func = zeros(11,1);

for kk = 1 : length(vec_kappa)
    kappa = vec_kappa(kk); % Intensity of the kick

% Floquet evolution 
U = expm(-1i*kappa*Jz*Jz/(2*j))*expm(-1i*p*Jy);      
    
otocm = zeros(1,N-tran);
for kicks = 1 : N
                      
    if (kicks > tran)          
       op1 = (U')*O1*U;
       otocm(kicks) = otoc((init')*init,O2,op1);
       O1 = op1;
    end
    
end

func(kk) = sum(otocm)/length(otocm);

end

plot(vec_kappa,func);

