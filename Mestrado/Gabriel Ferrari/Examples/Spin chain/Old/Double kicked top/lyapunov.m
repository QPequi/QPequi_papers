% Computes the Lyapunov spectrum of a dynamical map.

clear all;

% Global variables
kicks = 1000;
dim = 6; % dimention of the map
kappa = 5; % intensity of the kick
gamma = 0.01; % intensity of the interaction

% Initial conditions
x01 = 1;
y01 = 0;
z01 = 0;
x02 = 1;
y02  = 0;
z02 = 0;

lambda = zeros(dim,1); % Initializations
Q = eye(dim);

% Iteractions of the map - lyapunov exponents
for jj = 1 : kicks
    
% Kicked top map
Delta12 = kappa*x01 + gamma*x02;
Delta21 = kappa*x02 + gamma*x01;

xf1 = z01*cos(Delta12) + y01*sin(Delta12);
yf1 = -z01*sin(Delta12) + y01*cos(Delta12);
zf1 = -x01;
xf2 = z02*cos(Delta21) + y02*sin(Delta21);
yf2 = -z02*sin(Delta21) + y02*cos(Delta21);
zf2 = -x02;

% Jacobian
J11 = -kappa*(z01*sin(Delta12)-y01*cos(Delta12));
J12 = sin(Delta12);
J13 = cos(Delta12);
J14 = -gamma*(z01*sin(Delta12)-y01*cos(Delta12));
J15 = 0;
J16 = 0;

J21 = -kappa*(z01*cos(Delta12)+y01*sin(Delta12));
J22 = cos(Delta12);
J23 = -sin(Delta12);
J24 = -gamma*(z01*cos(Delta12)+y01*sin(Delta12));
J25 = 0;
J26 = 0;

J31 = -1;
J32 = 0;
J33 = 0;
J34 = 0;
J35 = 0;
J36 = 0;

J41 = -gamma*(z02*sin(Delta21)-y02*cos(Delta21));
J42 = 0;
J43 = 0;
J44 = -kappa*(z02*sin(Delta21)-y02*cos(Delta21));
J45 = sin(Delta21);
J46 = cos(Delta21);

J51 = -gamma*(z02*cos(Delta21)+y02*sin(Delta21));
J52 = 0;
J53 = 0;
J54 = -kappa*(z02*cos(Delta21)+y02*sin(Delta21));
J55 = cos(Delta21);
J56 = -sin(Delta21);

J61 = 0;
J62 = 0;
J63 = 0;
J64 = -1;
J65 = 0;
J66 = 0;

J = [J11 J12 J13 J14 J15 J16;
     J21 J22 J23 J24 J25 J26;
     J31 J32 J33 J34 J35 J36;
     J41 J42 J43 J44 J45 J46;
     J51 J52 J53 J54 J55 J56;
     J61 J62 J63 J64 J65 J66];
 
B = J*Q;

[Q,R] = qr(B);

lambda = lambda + log(diag(abs(R)));

x01 = xf1;
y01 = yf1;
z01 = zf1;
x02 = xf2;
y02 = yf2;
z02 = zf2;

end
lambda = lambda./kicks;
vec_max = max(lambda);