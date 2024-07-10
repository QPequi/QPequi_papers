%-------------------------------------------------------------------
% Initial parameters and definitions
l = 3;           % Angular momentum
dim = 2*l + 1;   % Dimension of the Hilbert space

% Angular momentum matrices in the z-representation (energy)
[Jx,Jy,Jz] = pauli(l);

Jp = Jx + 1i*Jy;
Jm = Jy - 1i*Jy;

% Energy eigenbasis and eigenvalues
eigenvec = flip(eye(dim)); % k-th eigenvector is eigenvc(:,k)
eigenval = abs(eig(Jz));   % k-th eigenvalue is engenval(k)

%-------------------------------------------------------------------
% Initial state
% Pure
%vec0 = (eigenvec(:,1) + eigenvec(:,dim))/sqrt(2);
%rho0 = vec0*(vec0');
% Thermal
beta = 0;
rho0 = diag(exp(-beta*eigenval));

%-------------------------------------------------------------------
% Initial entropy and coherence
sd0 = 0; % Diagonal entropy
for kk = 1 : dim
    if (rho0(kk,kk) < 0.0000001)
        sd0 = sd0 + 0;
    else    
    sd0 = sd0 - rho0(kk,kk)*log2(rho0(kk,kk));
    end
end
sd0 = real(sd0);

pop = eig(rho0);
svn0 = 0; % von Neumann entropy
for kk = 1 : dim
    if (pop(kk) < 0.0000001)
        svn0 = svn0 + 0;
    else    
    svn0 = svn0 - pop(kk)*log2(pop(kk));
    end
end
svn0 = real(svn0);

cor0 = sd0 - svn0; % Relative entropy of coherence

%-------------------------------------------------------------------
% Process
U = expm(1i*0.5*Jp + 1i*0.5*Jm);
rhof = U*rho0*(U');

%-------------------------------------------------------------------
% Final entropy and coherence
sdf = 0; % Diagonal entropy
for kk = 1 : dim
    if (rhof(kk,kk) < 0.0000001)
        sdf = sdf + 0;
    else    
    sdf = sdf - rhof(kk,kk)*log2(rhof(kk,kk));
    end
end
sdf = real(sdf);

pop = eig(rhof);
svnf = 0; % von Neumann entropy
for kk = 1 : dim
    if (pop(kk) < 0.0000001)
        svnf = svnf + 0;
    else    
    svnf = svnf - pop(kk)*log2(pop(kk));
    end
end
svnf = real(svnf);

corf = sdf - svnf; % Relative entropy of coherence
