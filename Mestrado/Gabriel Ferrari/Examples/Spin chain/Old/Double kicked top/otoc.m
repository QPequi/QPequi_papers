function f = otoc(state,V,W)
% OTOC for operators V and W (unitary ro Hermitian) with respect to state.
% W is the evolved operator, i.e. W = U'*W(0)*U

com = W*V - V*W;

f = trace(state*(com')*com);          
   
dig = 4;
f = round(f.*(10^dig))./(10^dig);
