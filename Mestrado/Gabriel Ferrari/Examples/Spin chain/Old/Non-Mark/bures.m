function f = bures(state,op)
% Bures angle between the density operators state and op.

aux = sqrtm(sqrtm(state)*op*sqrtm(state));

f = acos(trace(aux)*trace(aux));          
   
dig = 6;
f = round(f.*(10^dig))./(10^dig);
