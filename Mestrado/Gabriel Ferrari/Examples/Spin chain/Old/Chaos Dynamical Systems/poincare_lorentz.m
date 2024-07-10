% Integration of the equations for the duffung oscillator

global sigma r beta;

sigma = 10;
r = 28;
beta = 8/3;

step = 0.01;
time = 0 : step : 50;
init = [0 1 0];

% integration of the equations of motion
[t, xx] = ode45(@lorenz,time,init);

x = xx(:,1);
y = xx(:,2);
z = xx(:,3);

% time series
%{
figure(1)
subplot(2,1,1)
plot(t,x,'r')
subplot(2,1,2)
plot(t,p,'r')
ax = gca;
ax.FontSize = 22;
%}

% phase-space
figure(2)
plot3(x,y,z,'b')
ax = gca;
ax.FontSize = 22;

% Poincare section
nn = 1;
psi0 = 2*pi;
clear x1 x2;
for ii = 1 : length(x)
    
    pos1 = nn*psi0 - step/2;
    pos2 = nn*psi0 + step/2;
    if ( (t(ii) >= pos1) && (t(ii) <= pos2) )
       x1(nn)=x(ii);
       x2(nn)=y(ii);
       x3(nn)=z(ii);
       nn = nn + 1;
    end
    
end

figure(3)
plot3(x1,x2,x3,'b.')
ax = gca;
ax.FontSize = 22;