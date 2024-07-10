function op = tq_op(op1,op2,t1,t2,numq)
%--------------------------------------------------------------------------
% Build a two-qubit controlled gate in the form 
%
%     op(c,t) = I(2) X ... I(2) X op1_t1 X ... X op2_t2 X ... X I(2) 
%     t1 and t2 lables the sites of the two operators
%     I(2) is the identity operator in 2 dimensions
%     X means tensor product
% 
% The output is an operator over the entire Hilbert space of numq qubits
%-------------------------------------------

if nargin < 5
    error('Not enough input arguments');
end

if t1 == t2
    error('The two operators must act on different sites');
end

% Ordering t1 and t2
if t1 > t2
   t = t1;
   t1 = t2;
   t2 = t;
   aux = op1;
   op1 = op2;
   op2 = aux;
   clear t aux;
end

aux = cell(numq);

for ii = 1 : numq
    aux{ii} = eye(2);
    if (ii == t1)
        aux{ii} = op1;
    end
    if (ii == t2)
        aux{ii} = op2;
    end
end

op = kron(aux{1},aux{2});
for ii = 3 : numq
    op = kron(op,aux{ii});
end
        
end
