% Time vector
vect = 0:0.1:10;

% Angular momentun matrices in the Jz representation
[Jx,Jy,Jz] = pauli(1/2);

gamma0 = 0.1; % Intensity of Jx term.
% Initial Hamiltonian
H = Jz + gamma0*Jx;  
   
% Initial density matrix of the system
beta = 0.1;        % Inverse temperature (in units of the qubit gap) 
psi0 = expm(-beta*H)/trace(expm(-beta*H)); % Density mtrix (thermal state))

% time loop
for tt = 1 : length(vect) 

    gamma = gamma0*cos(vect(tt));
    % Hamiltonian
    H = Jz + gamma*Jx;

    % Evolved density matrix
    
    rhot = U*rho*U';
    
end
% Saving data
%save('SQ.dat','Entro','-ascii');

% Entropy production rate
% plot(vect(1:end-1),diff(Entro)./diff(vect));