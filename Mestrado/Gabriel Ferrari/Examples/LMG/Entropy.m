function s = Entropy(sigma)
% Compute the von Neumann entropy of the state sigma
%
%       S = Tr[ sigma log( sigma ) ]
%

[a,b] = size(sigma);
if a  ~= b
    disp('sigma is not a square matrix');
    s = 0;
    return;
end

D = eig(sigma);

s = 0;
for kk = 1 : a
    
    if D(kk) < 0.00000001
        s = s;
    else
        s = s - D(kk)*log(D(kk));
    end
    
end
s = real(s);