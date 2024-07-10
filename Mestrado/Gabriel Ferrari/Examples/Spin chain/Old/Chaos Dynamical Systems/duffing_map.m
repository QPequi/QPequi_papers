function xdot = duffing(t,x)
% forced duffing oscillator
%
%  d^2 x/dt^2 + gamma dx/dt + beta x^3 - x - (g/beta)cos(Omega t) = 0
%
%  The system can be rewritten in terms of two differential equations by
%  identifying p = dx/dt, which implies that d^2 xdt^2 = dp/dt
%
%  dx/dt = p
%  dp/dt = -gamma p - beta x^3 + x + (g/beta) cos(Ometa t) 

% gloval variables
global gamma Omega beta g

% system: p = x(2), x = x(1)
xdot(1) = x(2);
xdot(2) = -gamma*x(2) - beta*x(1)^3 + x(1) + (g/beta)*cos(Omega*t);

xdot=xdot';

end