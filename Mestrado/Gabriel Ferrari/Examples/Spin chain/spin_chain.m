%------------------------------------------
%
% Dynamics of spin chains
%
% Needs:
%       pauli.m
%       sq_gate.m
%       tq_gate.m
%------------------------------------------

[sx,sy,sz] = pauli(1/2);

N = 5;  % size of the chain
dim = power(2,N);

% Model
omega = 1.0;  % Bare frequency of the qubit
jxx = 1.0;  % Coupling constant in the xx direction
jyy = 1.0;  % Coupling constant in the yy direction
jzz = 0.5;  % Coupling constant in the zz direction

d = floor(N/2);

% Energy eigenvalues 
vect = 0 : 0.1 : 10;
epsd = 1.0;
alpha = 5;

% Interaction Hamiltonian
HI = zeros(dim,dim);
for ii = 1 : N-1
    HI = HI + (jxx/4)*tq_op(sx,sx,ii,ii+1,N) ...
            + (jyy/4)*tq_op(sy,sy,ii,ii+1,N) ...
            + (jzz/4)*tq_op(sz,sz,ii,ii+1,N);
end

eigenvalues = zeros(length(vect),dim);
for tt = 1 : length(vect)
    
    t = vect(tt);    
    eps = epsd*(tanh(alpha*(t-1)) + tanh(alpha));

% Free Hamiltonian
H0 = zeros(dim,dim);
for ii = 1 : N
    if (ii == d)
        H0 = H0 + ((omega+ eps)/2)*sq_op(sz,d,N);
    else
        H0 = H0 + (omega/2)*sq_op(sz,ii,N);
    end
end

H = H0 + HI;

[V,D] = eig(H);
[d,ind] = sort(diag(D));
eigenvalues(tt,:) = diag(D(ind,ind));
Vs = V(:,ind);
clear V D ind H0;
end

clear H;

 plot(t,eigenvalues(:,1))  

%psi0 = Vs(:,1); % Ground state

%------------------------------------------
% Dynamics
%------------------------------------------




%------------------------------------------
% Spectrum
%------------------------------------------
%{
% Ordered eigenvalues and eigenvectors
[V,D] = eig(H);
[d,ind] = sort(diag(D));
eigenvalues = diag(D(ind,ind));
Vs = V(:,ind);
clear V D ind;


%} 