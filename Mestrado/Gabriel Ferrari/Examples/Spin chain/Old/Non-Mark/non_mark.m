% Non-Markovianity measure for the kicked top
% The output of the program is the file non_mark.dat written in such a way
% that its collumns varies with kappa while its rows varies wiht j.

% Global parameters 
p = pi/2; % Energy of the kicked top
N = 1000; % Total number of kicks
tran = 100; % Transient to be disregarded

% Initial conditions uniformly distributed over a sphere
sample = 100; % Phase-space grid elements
[x,y,z] = unit(sample,1,1);
[phi,theta,radius] = cart2sph(x,y,z);
theta = sort(theta) + pi/2;
phi = sort(phi);

count = 0; % control vabriable

% System parameters
kappa = 0.0 : 0.5 : 5.0; % Intensity of the kick
ang = 10 : 10 : 100; % Angular momentum
kappa = kappa';
ang = ang';

% Data variable
mat_mark = zeros(length(ang),length(kappa));
mat_entro = zeros(length(ang),length(kappa));

teste = zeros(1,length(ang));

for jj = 1 : length(ang) 
 tic;   
j = ang(jj);
disp('j =  '); fprintf('%4d\n',j);    
    
dim = 2*j + 1; % Dimension of Hilbert space
  
[Jx,Jy,Jz] = pauli(j); % Angular momentum matrices in z basis
[sx,sy,sz] = pauli(1/2); % Pauli matrices

measure = zeros(1,length(kappa));
aux_entro = zeros(1,length(kappa));

for ii = 1 : length(kappa)
    
    measure_s = 0.0;
    entro_s = 0.0;
    
    for in = 1 : sample
        
    %  Orthogonal initial states
    th = theta(in);
    ph = phi(in);
    init0 = zeros(dim,1);
    init0(1) = 1;
    R = expm(-1i*th*(Jx*sin(ph) - Jy*cos(ph)));
    init = R*init0;
    
    [xn,yn,zn] = sph2cart(ph,th-pi/2,radius);
    [ph,th,radius] = cart2sph(-xn,-yn,-zn);
    th = th + pi/2;
    init_ort0 = zeros(dim,1);
    init_ort0(dim) = 1;
    R = expm(-1i*th(1)*(Jx*sin(ph(1)) - Jy*cos(ph(1))));
    init_ort = R*init0;

    teste = init'*init_ort;
    if teste > 0.00000001
        count = count + 1;
    end
    
    % Floquet evolution operator 
    U = expm(-1i*kappa(ii)*Jz*Jz/(2*j))*expm(-1i*p*Jy); 
  
    % auxiliary variables
    vecU1 = zeros(dim,1);
    vecU2 = zeros(dim,1);
    aux = 0.0; 
    
    vec1 = init;
    vec2 = init_ort;
    
    aux_entro_s = 0;
    
    % Time evolution
    for kicks = 1 : N
         
        % Transient
        if (kicks <= tran)
           vecU1 = U*vec1;
           clear vec1;
           vec1 = vecU1;
           clear vecU1;
           
           vecU2 = U*vec2;
           clear vec2;
           vec2 = vecU2;
           clear vecU2
        end
            
        if (kicks > tran)
            
           % Reduced states of the qubits
           av_sigmax = trace((vec1')*Jx*vec1)/j;  
           av_sigmay = trace((vec1')*Jy*vec1)/j;
           av_sigmaz = trace((vec1')*Jz*vec1)/j;
           rho_s1 = 0.5*(1 + av_sigmax.*sx + av_sigmay.*sy + av_sigmaz.*sz);
           
           av_sigmax = trace((vec2')*Jx*vec2)/j;  
           av_sigmay = trace((vec2')*Jy*vec2)/j;
           av_sigmaz = trace((vec2')*Jz*vec2)/j;
           rho_s2 = 0.5*(1 + av_sigmax.*sx + av_sigmay.*sy + av_sigmaz.*sz);
           
           % Norm before evolution
           op = rho_s1 - rho_s2;
           aux_normA = real(trace(sqrt(op*(op'))))/2; 
           clear rho_s1 rho_s2 op;
           
           % Time evolution
           vecU1 = U*vec1;
           clear vec1;
           vec1 = vecU1;
           clear vecU1;
           vecU2 = U*vec2;
           clear vec2;
           vec2 = vecU2;
           clear vecU2;
           
           % Reduced states of the qubits
           av_sigmax = trace((vec1')*Jx*vec1)/j;  
           av_sigmay = trace((vec1')*Jy*vec1)/j;
           av_sigmaz = trace((vec1')*Jz*vec1)/j;
           rho_s1 = 0.5*(1 + av_sigmax.*sx + av_sigmay.*sy + av_sigmaz.*sz);
           
           av_sigmax = trace((vec2')*Jx*vec2)/j;  
           av_sigmay = trace((vec2')*Jy*vec2)/j;
           av_sigmaz = trace((vec2')*Jz*vec2)/j;
           rho_s2 = 0.5*(1 + av_sigmax.*sx + av_sigmay.*sy + av_sigmaz.*sz);
           
           % Norm after evolution
           op = rho_s1 - rho_s2;
           aux_normB = real(trace(sqrt(op*(op'))))/2;
           
           % Considering only the positive values
           if (aux_normB - aux_normA > 0.0)
               aux = aux + aux_normB - aux_normA;
           end
           
           % Bures distance
           aux_entro_s = aux_entro_s + bures(rho_s1,eye(2));
     
           clear rho_s1 rho_s2 op;  
        end
            
    end % end of kicks loop
    
    measure_s = measure_s + real(aux);
    entro_s = entro_s + real(aux_entro_s/(N-tran));
    
    end   
    
    aux_entro(ii) = entro_s/sample;
    measure(ii) = measure_s/sample;
   
end % kappa loop

    if jj == 1
        mm = max(measure);
    end
    plot(kappa,measure./mm);
    hold on;
    
    mat_mark(jj,:) = measure;
    mat_entro(jj,:) = aux_entro;
    toc;
end % end j loop
hold off;

save('non_mark.dat', 'mat_mark','-ascii');
save('clausius.dat', 'mat_entro','-ascii');
save('kappa.dat', 'kappa','-ascii');
save('angular.dat', 'ang','-ascii');
init = [x,y,z];
save('data_init.dat','init','-ascii');