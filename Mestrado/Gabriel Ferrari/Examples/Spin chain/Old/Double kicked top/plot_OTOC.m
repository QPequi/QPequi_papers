clear all;

dataO = load('otoc_map.dat');
dataK = load('kappa.dat');
dataA = load('angular.dat');

m = max(max(dataO));

figure(1)
surf(dataA,dataK,dataO./m);
colormap(copper);
alpha(0.6);
shading interp;
axis([1 5 0 5 0 1]);
ax = gca;
ax.FontSize = 22;

figure(2)
pcolor(dataA,dataK,dataO./m);
colormap(copper);
shading interp;
axis([1 5 0 5]);
ax = gca;
ax.FontSize = 22;
