function [J1, J2, J3] = pauli(j)
%
%  Computes the angular momentum matrices in the z basis.
%  j is the total angular momentum
%

dim = 2*j+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Representation in Z basis
J1 = zeros(dim,dim);
J2 = zeros(dim,dim);
J3 = zeros(dim,dim);

% J1
aux1 = 0;
aux2 = 0;
for m1 = -j : j
    for m2 = -j : j
        
        if m1 == m2 + 1
            aux1 = sqrt(j*(j+1) - m2*(m2+1));
        end
        if m1 == m2 - 1
            aux2 = sqrt(j*(j+1) - m2*(m2-1));
        end

        pos1 = m1 + j + 1;
        pos2 = m2 + j + 1;
        
        J1(pos1,pos2) = (aux1 + aux2)/2.0;
        
        aux1 = 0;
        aux2 = 0;
        
    end
end

% J2
aux1 = 0;
aux2 = 0;
for m1 = -j : j
    for m2 = -j : j
        
        if m1 == m2 + 1
            aux1 = sqrt(j*(j+1) - m2*(m2+1));
        end
        if m1 == m2 - 1
            aux2 = sqrt(j*(j+1) - m2*(m2-1));
        end

        pos1 = m1 + j + 1;
        pos2 = m2 + j + 1;
        
        J2(pos1,pos2) = 1i*(aux1 - aux2)/2.0;
        
        aux1 = 0;
        aux2 = 0;
        
    end
end

% J3
m = j;
for count = 1 : dim
    J3(count,count) = m;
    m = m-1;
end