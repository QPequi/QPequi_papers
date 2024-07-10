% concurrence    Compute the concurrence of a two-qubit density matrix.

function c=concurrence(rho);
  y=[0 -i; +i 0];          
  R=(rho^0.5*kron(y,y)*conj(rho)*kron(y,y)*rho^0.5)^0.5;
  % Real is needed since MATLAB always adds a small imaginary part ...
  e=real(eig(R));
  e=-sort(-e);
  c=max(e(1)-e(2)-e(3)-e(4),0);