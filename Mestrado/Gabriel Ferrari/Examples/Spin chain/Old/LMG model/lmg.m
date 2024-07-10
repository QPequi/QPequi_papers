% Global parameters
j = 500;          % Angular momentum
dim = 2*j + 1;    % Dimension of the system
gammax = 1;       % Jx^2 coupling 
gammay = 0;       % Jy^2 coupling

% Angular momentun matrices in the Jz representation
[Jx,Jy,Jz] = pauli(j);

%-------------------------------------------------------------
%
%      LEVEL SPACE STATISTICS
%
%-------------------------------------------------------------
%{
vech = 0:0.2:2;

for hh = 1 : length(vech)
    
h = vech(hh);
H = -h*Jz - (gammax*mpower(Jx,2) + gammay*mpower(Jy,2))/(2*j);

% Diagonalization
[V,D] = eig(H);
[d,ind] = sort(diag(D));
Ds = D(ind,ind);
Vs = V(:,ind);
vec_eig = zeros(dim,1);
for kk = 1 : dim
    vec_eig(kk) = Ds(kk,kk);
end

% Testing unitarity
teste = max(max(abs(Vs*(Vs') - eye(dim))));
if (teste > 1e-12) 
    disp('error in diagonalization'); 
end

clear V D d ind H;

% Normalizing the eigenvalues
vec_eig = vec_eig + abs(vec_eig(1));
vec_eig = vec_eig/vec_eig(dim);

% level spacing statistics
spacing = zeros(dim-1,1);
for kk = 1 : dim-1
    spacing(kk) = vec_eig(kk+1) - vec_eig(kk);
end

bin = 100;
num_bins = sort(linspace(1,0,bin));
dist = zeros(length(num_bins),length(vech));
for ii = 1 : length(num_bins)-1
    count = 0;
    for kk = 1 : dim
        a = vec_eig(kk);
        if ( (a >= num_bins(ii)) && (a < num_bins(ii+1)))
            count = count + 1; 
        end
    end    
    dist(ii,hh) = count;
end
dist(end,:) = [];

end
dist = sort(dist,'descend')/max(max(dist));
save('level.dat','dist','-ascii')
% Each column contains the statistics associated with some value of h
%}

%-------------------------------------------------------------
%
%      INVARIANT ENTROPY
%
%-------------------------------------------------------------
%{
h = 0;
H0 = -h*Jz - (gammax*mpower(Jx,2) + gammay*mpower(Jy,2))/(2*j);

% Diagonalization and initial state
[V,D] = eig(H0);
[d,ind] = sort(diag(D));
Ds = D(ind,ind);
Vs = V(:,ind);
psi0 = Vs(:,1); % Ground state

% Initial density matrix of the system 
%beta = 0.1; % In units of gammax
%psi0 = expm(-beta*H0)/trace(expm(-beta*H0)); % Thermal density mtrix

vech = [0 0.2 0.5 0.8 1.2];
vect = 0:0.1:50;
entro = zeros(length(vect),length(vech));
entro_av = zeros(length(vech),1);
entro_diff = zeros(length(vect)-1,length(vech));
for hh = 1 : length(vech)

h = vech(hh);

HF = -h*Jz - (gammax*mpower(Jx,2) + gammay*mpower(Jy,2))/(2*j);

aux = zeros(length(vect),1);

for tt = 1 : length(vect)
    
     % Evolved state
     U = expm(-1i*HF*vect(tt));
     psit = U*psi0;
     rhot = psit*(psit');
     
     % Changing basis of the density matrix
     rho = (Vs')*rhot*(Vs);
     
     % Entropy
     s = 0.0;
     for jj = 1 : dim
         if (rho(jj,jj) < 1e-6)
             s = s + 0;
         else
             s = s - rho(jj,jj)*log(rho(jj,jj));
         end
     end
     
     aux(tt) = real(s);
     
end 
entro(:,hh) = aux; % Entropy
aux1 = diff(aux);
entro_diff(:,hh) = aux1; % Entropy rate
entro_av(hh) = trapz(aux)/vect(length(vect)); % Average entropy

end
save('time.dat','vect','-ascii')
save('entro.dat','entro','-ascii')
save('entro_av.dat','entro_av','-ascii')
save('entro_diff.dat','entro_diff','-ascii')

%}

%-------------------------------------------------------------
%
%       LOSCHMIDT ECHO
%
%-------------------------------------------------------------
%{
h = 0;
H0 = -h*Jz - (gammax*mpower(Jx,2) + gammay*mpower(Jy,2))/(2*j);

% Diagonalization and initial state
[V,D] = eig(H0);
[d,ind] = sort(diag(D));
Ds = D(ind,ind);
Vs = V(:,ind);
psi0 = Vs(:,1); % Ground state

% Initial density matrix of the system 
%beta = 0.1; % In units of gammax
%psi0 = expm(-beta*H0)/trace(expm(-beta*H0)); % Thermal density mtrix

vech = [0 0.2 0.5 0.8 1.2];
vect = 0:0.1:50;
echo = zeros(length(vect),length(vech));
for hh = 1 : length(vech)
tic;
     h = vech(hh);

     HF = -h*Jz - (gammax*mpower(Jx,2) + gammay*mpower(Jy,2))/(2*j);

     % Evolved state
     for tt = 1 : length(vect)
        U = expm(-1i*HF*vect(tt));
        psit = U*psi0;
     
        echo(tt,hh) = power(abs((psit')*psi0),2);
     end
toc;
end

save('echo.dat','echo','-ascii')
%}