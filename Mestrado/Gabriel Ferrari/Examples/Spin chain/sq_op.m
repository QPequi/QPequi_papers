function op = sq_op(op1,t,numq)
%--------------------------------------------------------------------------
% Build a single qubit gate in the form 
%
%     op = I(2) X ... I(2) X ... X op1_t X ... X I(2) 
%     t indicates the site where the operator op1 acts
%     I(2) is the identity operator in 2 dimensions
%     X means tensor product
%
% The output is an operator over the entire Hilbert space of numq qubits
%--------------------------------------------------------------------------

if nargin < 3
    error('Not enough input arguments');
end

if (t == 1)
        
    aux = op1;
        
    for ll = 2 : numq      
        aux = kron(aux,eye(2));
    end
        
end   
    
if (t > 1)
        
   aux = eye(2);
        
   for ll = 1 : (t - 2)      
       aux = kron(eye(2),aux);
   end
    
   aux = kron(aux,op1);
        
   for ll = (t + 1) : numq
       aux = kron(aux,eye(2));
   end
        
end

op = aux;

end