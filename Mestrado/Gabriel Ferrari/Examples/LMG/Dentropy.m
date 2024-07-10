function s = Dentropy(V,sigma)
% Compute the diagonal entropy of the state sigma, which is the Shannon
% entropy of the probability distribution contained in the diagonal of the
% density matrix. Diagonal with respect to the basis defined in V
%
%       rho = (V')*sigma*V
%       S = Sum_k rho(kk,kk) log( rho(kk,kk) )
%
% It is assumed that the state sigma is written in the energy basis

[a,b] = size(sigma);
if a  ~= b
    disp('sigma is not a square matrix');
    s = 0;
    return;
end

rho = (V')*sigma*V;

s = 0;
for kk = 1 : a
    
    if rho(kk,kk) < 0.00000001
        s = s;
    else
        s = s - rho(kk,kk)*log(rho(kk,kk));
    end
    
end
s = real(s);
        