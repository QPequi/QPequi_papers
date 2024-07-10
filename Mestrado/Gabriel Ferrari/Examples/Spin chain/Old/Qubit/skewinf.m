% skewinf   Wigner-Yanase skew information. Usage: skewinf(rho,A), where
%           rho is a density matrix and A is an operator.

function s=skewinf(rho,A)

s=trace(A^2*rho)-trace(rho^0.5*A*rho^0.5*A);