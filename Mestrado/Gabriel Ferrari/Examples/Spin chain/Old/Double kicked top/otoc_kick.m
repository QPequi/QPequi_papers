%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The double kicked top. Indexes 1 and 2 characterize each indivisual
% system. All the matrices are writen in the z basis of each system.

clear all;

% Global parameters 

% Energy of the kicked top
p1 = pi/2; 
p2 = pi/2;

% Number of kicks
N1 = 1000;
N2 = 1000;

 % Pauli matrices
[sigmax,sigmay,sigmaz] = pauli(1/2);

% Angular momentum and dimension of the Hilbert space
j1 = 7/2;
j2 = 7/2;
dim1 = 2*j1 + 1;
dim2 = 2*j2 + 1;
  
% Angular momentum matrices
[Jx1,Jy1,Jz1] = pauli(j1);
[Jx2,Jy2,Jz2] = pauli(j2);

% Global initial state vector - thermal
beta = 0.1;

init1 = zeros(dim1,dim1);
kk = 1;
z = 0.0;
for m = -j1 : 1 : j1
    init1(kk,kk) = exp(-beta*m);
    kk = kk + 1;
    z = z + exp(-beta*m); 
end
init1 = init1./z;

init2 = zeros(dim2,dim2);
kk = 1;
z = 0.0;
for m = -j2 : 1 : j2
    init2(kk,kk) = exp(-beta*m);
    kk = kk + 1;
    z = z + exp(-beta*m); 
end
init2 = init2./z;

init = kron(init1,init2);

clear init1 init2 z kk;

% OTOC operators
V = expm(-1i*p1*kron(Jy1,eye(dim2)));
W = expm(-1i*p1*kron(Jy1,eye(dim2))); 

% Interaction between both systems
gamma = 0.00;

% Range for the intensity of the kicks. They can be changed individually
% latter.
vec_kappa = 0 : 1 : 6;

data_otoc = zeros(length(vec_kappa),1);
data_ent = zeros(length(vec_kappa),1);

for kk = 1 : length(vec_kappa)
    
    kappa1 = vec_kappa(kk); % Intensity of the kick - first system
    kappa2 = vec_kappa(kk); % Intensity of the kick - second system

% Total Floquet evolution
u1 = expm(-1i*kappa1*Jz1*Jz1/(2*j1))*expm(-1i*p1*Jy1);
u2 = expm(-1i*kappa2*Jz2*Jz2/(2*j2))*expm(-1i*p2*Jy2);
ui = expm(-1i*gamma*kron(Jz1,Jz2)/j1);
U = ui*kron(eye(dim1),u2)*kron(u1,eye(dim2));

% Auxiliary variables
otocm = zeros(N1,1);
rho0 = init;

% Main loop
for kicks = 1 : N1   
      % Evolved state
      rhou = U*rho0*(U'); 
      rho0 = rhou;
      
       % OTOC
       Wu = (U')*W*U;
       otocm(kicks) = otoc(init,V,Wu);
       W = Wu;
end

% Average OTOC
data_otoc(kk) = sum(otocm)/length(otocm);

% Entanglement
rho_red1 = ptrace(init,2,[dim1,dim2]);
rho_red2 = ptrace(init,1,[dim1,dim2]);
data_ent(kk) = entropy(rho_red1) + entropy(rho_red2) - entropy(init);

end

figure(1)
plot(vec_kappa,data_otoc);

figure(2)
plot(vec_kappa,data_ent);