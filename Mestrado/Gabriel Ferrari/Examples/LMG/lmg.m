clear all

for j = 10 : 10 : 50
    disp(j)
dim = 2*j + 1;          % dimension of the Hilbert space
[jx,jy,jz] = pauli(j);  % angular momentum matrices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spectrum of the Hamiltonian
%{
lvec = 0 : 0.01: 0.5;

energies = zeros(length(lvec),dim);

for ii = 1 : length(lvec)
    lambda = lvec(ii);
    H0 = (lambda/j)*mpower(jx,2) - jz; 
    spec = eig(H0);
    energies(ii,:) = spec;
end

for ii = 1 : dim
plot(lvec,energies(:,ii))
hold on
end
hold off
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time evolution

% Initial state - ground state of the initial Hamiltonian
% Since lambda(0) = 0
H0 = -jz;
[V,D] = eig(H0); % eigenvector k = V(:,k)
state0 = V(:,1); % ground state is the initial state of the system
clear V D;

alpha = 5;
tspan = [0 3];   % Integration time
y0 = state0;     % Initial conditions for the integration

%opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,y] = ode45(@(t,y) sch(t,y,alpha,j),tspan,y0); % Integration
% The k-th solution - coefficient of the expansion - are given by y(:,k)

%plot(t,abs(y(:,3)))

% density operator
coeft = zeros(length(t),dim);
rhot = cell(length(t),1);

[V,D] = eig(jz); % H0 = Jz
for ii = 1 : length(t)
    state = zeros(dim,1);
    for jj = 1 : dim
        state = state + y(ii,jj)*V(:,jj);
    end
    % rhot{ii} contains the density operator at time t = t(ii);
    rhot{ii} = state*(state');
end

% Testing if the trace of rho is one.
%{
teste = zeros(length(t),1);
for tt = 1 : length(t)
    teste(tt) = real(trace(rhot{tt}));
end
plot(t,teste)
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Total energy variation of the system
aux = zeros(length(t),1);
for kk = 1 : length(t)
    aux(kk) = trace(rhot{kk}*mpower(jx,2));
end
du = (alpha/j)*(sech(alpha*(t-1)).*sech(alpha*(t-1))).*aux;

% Work
lambda = (1/(4*j))*(tanh(alpha*(t-1)) + tanh(alpha)); % Protocol
energies = zeros(length(t),dim);
for kk = 1 : length(t) 
    H = lambda(kk)*mpower(jx,2) - jz;  % Hamiltonian at time t
    [V,D] = eig(H);                    % Diagonalization 
    energies(kk,:) = diag(D);          % Energies
end

denergies = zeros(length(t),dim);
for kk = 1 : dim 
    denergies(:,kk) = gradient(energies(:,kk))./gradient(t);
end

dwork = zeros(length(t),1);
basis = cell(length(t),1); %energy basis
for kk = 1 : length(t)
    H = lambda(kk)*mpower(jx,2) - jz;  % Hamiltonian at time t
    [V,D] = eig(H);  % Diagonalization 
    
    basis{kk} = V;  % Keep the energy eigenbasis for every t
    
    hdot = diag(denergies(kk,:));
    
    dwork(kk) = trace(rhot{kk}*V*hdot*(V'));
end

dheat = du - dwork;   % heat variation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total energies

u = zeros(length(t),1);
work = zeros(length(t),1);
heat = zeros(length(t),1);
clear vec dt;
for kk = 1 : length(t)
    % Differentials of time
    dt(1)=0;
    for nn = 2: kk
       dt(nn) = t(nn) - t(nn-1);
    end
    
    vec = du(1:kk).*(dt');
    u(kk) = real(trapz(vec)); % Energy
    
    vec = dwork(1:kk).*(dt');
    work(kk) = real(trapz(vec)); % Work
    
    vec = dheat(1:kk).*(dt');
    heat(kk) = real(trapz(vec)); % Heat
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coherences
coh = zeros(length(t),1);
for kk = 1 : length(t)
    coh(kk) = Dentropy(basis{kk},rhot{kk}) - Entropy(rhot{kk});
end

plot(t,coh,'k','LineWidth',2);
%hold on
%plot(t,work,'--r','LineWidth',2);
%hold on
%plot(t,heat,'.b','LineWidth',2);
hold on

end
hold off
xticks([0 1 2 3 4 5])
%xticklabels({'0','$\pi$','$2\pi$'})
xlabel('$t$','Interpreter','latex')
ylabel('$W_{\mathrm{inv}}$','Interpreter','latex')
set(gca,'FontSize',33)