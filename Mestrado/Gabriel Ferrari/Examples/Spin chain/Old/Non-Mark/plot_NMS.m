clear all;

dataNM = load('non_mark.dat');
dataS = load('clausius.dat');
dataK = load('kappa.dat');

m = max(dataNM(1,:));

figure(1)
plot(dataK,dataNM(1,:)./m,'+', 'MarkerSize', 8)
hold on
plot(dataK,dataNM(2,:)./m, 'o', 'MarkerSize', 8)
hold on
plot(dataK,dataNM(3,:)./m, '*', 'MarkerSize', 8)
hold on
plot(dataK,dataNM(4,:)./m, 'x', 'MarkerSize', 8)
hold on
plot(dataK,dataNM(5,:)./m, 's', 'MarkerSize', 8)
hold on
plot(dataK,dataNM(6,:)./m, 'd', 'MarkerSize', 8)
hold on
plot(dataK,dataNM(7,:)./m, 'p', 'MarkerSize', 8)
hold on
plot(dataK,dataNM(8,:)./m, 'h', 'MarkerSize', 8)
hold on
plot(dataK,dataNM(9,:)./m, '^', 'MarkerSize', 8)
hold on
plot(dataK,dataNM(10,:)./m, 'v', 'MarkerSize', 8)
hold off
axis([0 5 0 1]);
ax = gca;
ax.FontSize = 22;

m = max(dataS(1,:));

figure(2)
plot(dataK,dataS(1,:)./m, '+', 'MarkerSize', 8)
hold on
plot(dataK,dataS(2,:)./m, 'o', 'MarkerSize', 8)
hold on
plot(dataK,dataS(3,:)./m, '*', 'MarkerSize', 8)
hold on
plot(dataK,dataS(4,:)./m, 'x', 'MarkerSize', 8)
hold on
plot(dataK,dataS(5,:)./m, 's', 'MarkerSize', 8)
hold on
plot(dataK,dataS(6,:)./m, 'd', 'MarkerSize', 8)
hold on
plot(dataK,dataS(7,:)./m, 'p', 'MarkerSize', 8)
hold on
plot(dataK,dataS(8,:)./m, 'h', 'MarkerSize', 8)
hold on
plot(dataK,dataS(9,:)./m, '^', 'MarkerSize', 8)
hold on
plot(dataK,dataS(10,:)./m, 'v', 'MarkerSize', 8)
hold off
axis([0 5 0.2 1]);
ax = gca;
ax.FontSize = 22;
