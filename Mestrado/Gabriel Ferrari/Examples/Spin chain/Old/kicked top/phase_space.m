% This program computes the phase-space portraint for the kicked top.

clear all;

% Initial conditions uniformly distributed over a sphere
sample = 300; % Phase-space grid elements
[x,y,z] = unit(sample,1,1);
[phi,theta,radius] = cart2sph(x,y,z);
theta = -sort(theta);
phi = sort(phi);

kappa = 1;  % Intensity of the kicks
kicks = 10000; % Number of kicks  

[x,y,z] = sph2cart(phi,-theta,radius);

xf = zeros(kicks,1);
yf = zeros(kicks,1);
zf = zeros(kicks,1);

phif = zeros(kicks,1);
thetaf = zeros(kicks,1);

figure(1);
for ii = 1 : sample

    x0 = x(ii);
    y0 = y(ii);
    z0 = z(ii);

for ll = 1 : kicks
% Kicked top map
xf(ll) = z0*cos(kappa*x0) + y0*sin(kappa*x0);
yf(ll) = -z0*sin(kappa*x0) + y0*cos(kappa*x0);
zf(ll) = -x0;

x0 = xf(ll);
y0 = yf(ll);
z0 = zf(ll);
end

% 3D plot
subplot(2,1,1)
plot3(xf,yf,zf,'.','MarkerSize',2,'MarkerEdgeColor','blue',...
        'MarkerFaceColor','blue');
hold on

% 2D plot
[phif,thetaf,radf] = cart2sph(xf,yf,zf);
thetaf = thetaf + pi/2;
subplot(2,1,2)
plot(phif,thetaf,'.','MarkerSize',2);
hold on

end
plot(-pi/2,pi/2,'o','MarkerSize',15);
hold on;
axis([-pi pi 0 pi]);

subplot(2,1,1)
[xs,ys,zs] = sphere(100);
fig = surf(xs,ys,zs,'FaceAlpha',0.1,'EdgeColor', 'none');
set(fig,'FaceColor', 'yellow');
hold off