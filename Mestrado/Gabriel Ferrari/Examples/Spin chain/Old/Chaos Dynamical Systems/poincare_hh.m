function g=henon_heiles(energy,tmax,niter)
%----------------------------------------------------------------
% Calculates the orbits of the Henon-Heiles Hamiltonian and the 
% intersections with a Poincare section
% 
%  Call henon_heiles(energy,tmax,niter)
%
%  Input parameters
%     0 < energy < 1/6
%     tmax = integration time
%     niter = integer > 1 = no of different orbits
% 
%  Example: henon_heiles(1/10, 500,15)
%  
% The Poincare section is at x = 0.
% The phase space plot in variables (y,py)
% For given (E, y, py) there are two choices for px.
% ----------------------------------------------------------------

close all,
E = energy;
timespan = [0  tmax];
e1 = E*6;

if e1 > 1
disp('>> Energy is too high, motion unbounded!');
return
end 
 
s_1 = linspace(eps, 1-eps, niter);
p_1 = sqrt(2*E)*sin(s_1*pi/2);
p_2 = sqrt(2*E)*cos(s_1*pi/2);

zz=[];
%%----------------------------------------------------------------
for iter=1:niter

iv = [0, p_1(iter),0 , p_2(iter)]';
 
options=odeset('AbsTol',1e-10,'RelTol',1e-5,'Events',@events );
[T,Y, TE,YE,IE]=ode45(@ff,timespan,iv,options);

zz=[zz; YE(:,3:4)];
disp(iter),

end 
%%----------------------------------------------------------------

figure(1)
plot(zz(:,1), zz(:,2),'.');
hold on;
plot(zz(:,1), - zz(:,2),'.');
hold off;

%% ---------------------------------------------------------------
function ydot=ff(t,y)
ydot = [y(2); 
	   - y(1)*(1+ 2*y(3));
	   y(4);
        - y(1)^2 - y(3) + y(3)^2];
% %%--------------------------------------------------------------
% function g=energy(y);
% g=0.5*sum(y.^2) + y(3)*(y(1)^2 - (y(3)^2)/3);

% %%--------------------------------------------------------------
function  [value,isterminal,direction] =  events(t,y)
global rho
value = y(1);
isterminal = 0;
direction = 1;
%%----------------------------------------------------------------
