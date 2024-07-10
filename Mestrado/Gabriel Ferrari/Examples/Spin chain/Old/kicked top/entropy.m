function f = entropy(state)
% computes the diagonal entropy of a system
% state must the de density operator in the energy eigenbasis

n = length(state);

f = 0;
for ii = 1 : n
    if state(ii) ~= 0
        f = f - state(ii)*log2(state(ii));
    else
        f = f + 0;
    end
end

end

