% Computation of the averaged Husimi distribution. Average over the kicks
% given a specific initial condition. Also Wehrl entropy

% Global parameters 
p = pi/2; % Energy of the kicked top
N = 1000; % Total number of kicks
tran = 100; % Transient to be disregarded

% Initial conditions
sample = 100; % Phase-space grid elements
[x,y,z] = unit(sample,1,1);
[phi,theta,rad] = cart2sph(x,y,z);
theta = theta + pi/2;
phi = sort(phi);
theta = sort(theta);
step_theta = 0;
step_phi = 0;
for ltheta = 1 : sample - 1
    step_theta = step_theta + theta(ltheta+1) - theta(ltheta);
end    
for lphi = 1 : sample - 1 
    step_phi = step_phi + phi(lphi+1) - phi(lphi);
end 
step_theta = step_theta/sample;
step_phi = step_phi/sample;

%[sigmax,sigmay,sigmaz] = Pauli(1/2); % Pauli matrices in z basis
kappa = 0.0 : 0.05 : 5.0; % Intensity of the kick
j = 1 : 1 : 200; % Angular momemtum

fidw = fopen('wehrl.dat','w'); % Entropy data

for count_j = 1 : 1 : length(j) % Loop for the angular momentum
  
fprintf('j = %4.2f\n', j(count_j));
dim = 2*j(count_j) + 1; % Dimension of Hilbert space

[Jx,Jy,Jz] = pauli(j(count_j)); % Angular momentum matrices in z basis

% Initial density operator
init = zeros(dim,1);
init(1) = 1;

for count_k = 1 : 1 : length(kappa) % Loop for the intensity of the kick
    
    fprintf('     kappa = %8.2f\n', kappa(count_k));

    fidh = fopen('husimi.dat', 'w'); % Husimi data file

    % Floquet evolution and average Husimi distribution
    U = expm(-1i*kappa(count_k)*Jy*Jy/(2*j(count_j)))*expm(-1i*p*Jz);
    rhom = init;
    rhomu = zeros(dim,1);
    aux = zeros(dim,1);
        
    % time-averaged state
    for kicks = 1 : N
            
        if (kicks <= tran)
           rhomu = U*rhom;
           rhom = rhomu;
        end
            
        if (kicks > tran)
           rhomu = U*rhom;
           aux = aux + rhomu;
           rhom = rhomu;
        end
            
    end
    rhom = aux/(N-tran);    
        
    % Husimi function
    for ltheta = 1 : sample
       for lphi = 1 : sample
            
       % Spin coherent state representation
       R = expm(-1i*theta(ltheta)*(Jx*sin(phi(lphi)) - Jy*cos(phi(lphi))));
       rhor = R*init;
       dist = trace((rhom')*rhor); % Husimi distribution
       
       fprintf(fidh, '%10.6f     %10.6f     %10.6f\n',...
                phi(lphi), theta(ltheta), dist);  
        
        end % theta loop 
 
    end % phi loop
    fclose(fidh);

% Wehrl entropy
data = load('husimi.dat');

aux = zeros(length(data(:,1)),length(data(:,2)));
func = data(:,3);

for ii = 1 : length(data(:,1))
    kk = 1;
    for jj = 1 : length(data(:,2))
        
        if (func(kk)< 0.0001)
            aux(ii,jj) = 0.0;
        else
            aux(ii,jj) = -func(kk)*log(func(kk));
        end
        kk = kk + 1;
        
    end
end

% Integration
entro = ((2*j(count_j) + 1)/(4*pi))*trapz(trapz(aux))*step_theta*step_phi; % double integration

fprintf(fidw, '%3.1f     %4.1f     %8.4f\n', kappa(count_k),...
                j(count_j), entro); 

end % End of the loop for the intensity of the kick

end % End of the angular momentum loop

fclose(fidw);