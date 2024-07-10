% Computes the Lyapunov spectrum of a dynamical map.

kicks = 100;
dim = 3; % dimention of the map
kappa = 5; % intensity of the kick

%for ll =  1 : length(kappa)
    
init = [0 0 0]';    

lambda = zeros(dim,1); % Initializations
Q = eye(dim);
% Iteractions of the map - lyapunov exponents
for jj = 1 : kicks
    
% Kicked top map
xf = z0*cos(kappa(ll)*x0) + y0*sin(kappa(ll)*x0);
yf = -z0*sin(kappa(ll)*x0) + y0*cos(kappa(ll)*x0);
zf = -x0;

% Jacobian and QR decomposition at the evolved point
J = [-kappa(ll)*(z0*sin(kappa(ll)*x0)-y0*cos(kappa(ll)*x0)) sin(kappa(ll)*x0) cos(kappa(ll)*x0);...
      -kappa(ll)*(z0*cos(kappa(ll)*x0)+y0*sin(kappa(ll)*x0)) cos(kappa(ll)*x0) -sin(kappa(ll)*x0);...
                          -1                       0              0      ];
B = J*Q;

[Q,R] = qr(B);

lambda = lambda + log(diag(abs(R)));

x0 = xf;
y0 = yf;
z0 = zf;

end
lambda = lambda./kicks;
vec_max(ii) = max(lambda);

result = 0.0;
for kk = 1:dim
    if (lambda(kk) > 0)
        result = result + lambda(kk);
    end
end
vec_sum(ii) = result;

end % end initial conditions loop

lambda_sum(ll) = sum(vec_sum)/length(vec_sum); % Kolmogovo-Sinai entropy
lambda_max(ll) = sum(vec_max)/length(vec_max); % Maximum exponent

end % end kappa loop

figure(1)
plot(kappa,lambda_max,'.','MarkerSize',4);

fid = fopen('lyap_coef.dat','w');
for ii = 1 : length(lambda_sum)
    fprintf(fid,'%10.4f      %10.4f\n', kappa(ii),lambda_max(ii));
end
fclose(fid);