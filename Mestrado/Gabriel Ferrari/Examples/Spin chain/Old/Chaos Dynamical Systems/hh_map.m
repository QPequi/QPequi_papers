function f=hh_map(t,X)
%
% This function provides the dynamical equations for Henon-Heiles system,
% as well as the linearized systems and the variational equations for the
% computation of the Lyapunov exponents
%
%       dx/dt  = px
%       dy/dt  = py
%       dpx/dt = -x - 2*lambda*x*y
%       dpy/dt = -y -lambda*(x^2 - y^2)
%

% Values of parameters
global lambda;

x = X(1); 
y = X(2); 
px = X(3); 
py = X(4);

Y= [X(5), X(9),  X(13), X(17);
    X(6), X(10), X(14), X(18);
    X(7), X(11), X(15), X(19);
    X(8), X(12), X(16), X(20)];

f=zeros(20,1);

%Henon-Heiles equations
f(1) = px;
f(2) = py;
f(3) = -x - 2*lambda*x*y;
f(4) = -y - lambda*(x^2 - y^2);

%Linearized system
 Jac=[       0,            0,       1, 0;
             0,            0,       0, 1;
      -1-2*lambda*y,  -2*lambda*x,  0, 0;
        -2*lambda*x, -1+2*lambda*y, 0, 0];
  
%Variational equation   
f(5:20)=Jac*Y;

%Output data must be a column vector


