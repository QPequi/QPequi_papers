function xdot = lorenz(t,x)
%
%  Lorenz equation 
%
%               dx/dt = SIGMA*(y - x)
%               dy/dt = R*x - y -x*z
%               dz/dt= x*y - BETA*z
%
%        In demo run SIGMA = 10, R = 28, BETA = 8/3
%        Initial conditions: x(0) = 0, y(0) = 1, z(0) = 0;
%        Reference values for t=10 000 : 
%              L_1 = 0.9022, L_2 = 0.0003, LE3 = -14.5691
%
%        See:
%    K. Ramasubramanian, M.S. Sriram, "A comparative study of computation 
%    of Lyapunov spectra with different algorithms", Physica D 139 (2000) 72-86.
%
% --------------------------------------------------------------------
% Copyright (C) 2004, Govorukhin V.N.
% Values of parameters

global sigma r beta;

xdot(1) = sigma*(x(2) - x(1));
xdot(2) = -x(1)*x(3) + r*x(1) - x(2);
xdot(3) = x(1)*x(2) - beta*x(3);

xdot = xdot';
%Output data must be a column vector