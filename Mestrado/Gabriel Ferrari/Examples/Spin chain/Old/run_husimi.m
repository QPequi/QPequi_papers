% Global parameters
j = 100;            % Angular momentum
dim = 2*j + 1;    % Dimension of the system
gammax = 1;       % Jx^2 coupling 
gammay = 0;       % Jy^2 coupling 
vect = 0:0.1:10;  % Time vector

% Angular momentun matrices in the Jz representation
[Jx,Jy,Jz] = pauli(j);

h = 0; % Initial Hamiltonian
H0 = -h*Jz - (gammax*mpower(Jx,2) + gammay*mpower(Jy,2))/(2*j);

% Initial density matrix of the system 
beta = 0.1; % In units of gammax
psi0 = expm(-beta*H0)/trace(expm(-beta*H0)); % Thermal density mtrix
    
% Final Hamiltonian
h = 0.8;
HF = -h*Jz - (gammax*mpower(Jx,2) + gammay*mpower(Jy,2))/(2*j);

% Change of basis
[V,D] = eig(HF);
[d,ind] = sort(diag(D));
Ds = D(ind,ind);
Vs = V(:,ind);
clear V D d ind;

% Testing unitarity
teste = max(max(abs(Vs*(Vs') - eye(dim))));
if (teste > 1e-12) 
    disp('error in diagonalization'); 
end

entro = zeros(length(vect),1); % Entropy data vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   INVARIANT ENTROPY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for tt = 1 : length(vect)
    
     % Evolved state
     U = expm(-1i*HF*vect(tt));
     psit = U*psi0*(U');
     
     % Changing basis of the density matrix
     rho = (Vs')*psit*(Vs);
     
     % Entropy
     s = 0.0;
     for jj = 1 : dim
         if (rho(jj,jj) < 1e-6)
             s = s + 0;
         else
             s = s - rho(jj,jj)*log(rho(jj,jj));
         end
     end
     
     entro(tt) = real(s);
     
end     
plot(vect,entro)       
     % Integration over theta and phi 
 %    Entro(tt) = real(dim*trapz(theta,trapz(phi,sq))/(4*pi));   
%end

% Saving data
%save('SQ.dat','Entro','-ascii');

% Entropy production rate
% plot(vect(1:end-1),diff(Entro)./diff(vect));