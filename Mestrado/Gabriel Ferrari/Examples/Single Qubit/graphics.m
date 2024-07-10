% Readind data
coh_cos = load('coh_cos.txt');
heat_cos = load('heat_cos.txt');
work_cos = load('work_cos.txt');

coh_cos_the = load('coh_cos_the.txt');
heat_cos_the = load('heat_cos_the.txt');
work_cos_the = load('work_cos_the.txt');

entroC_cos_the = load('entroC_cos_the.txt');
entroF_cos_the = load('entroF_cos_the.txt');
entroS_cos_the = load('entroS_cos_the.txt');

% Time vector - Same for all
time = work_cos(:,1);

% Pure state |1>
figure(1) 
plot(time,coh_cos(:,2),'k');
hold on
plot(time,work_cos(:,2)./2,'--r');
hold on
plot(time,heat_cos(:,2),':b');
hold off
xlabel('$t$', 'FontSize', 50, 'Interpreter','latex');
ylabel('$W_{\mathrm{inv}},\, Q_{\mathrm{inv}},\, C$', 'FontSize', 50, 'Interpreter','latex');
set(gca,'FontSize', 38, 'TickLabelInterpreter','latex')
set(gca,'TickLabelInterpreter','latex')

set(gca,'FontSize', 38, 'TickLabelInterpreter','latex')
set(gca,'TickLabelInterpreter','latex')

% Thermal state - beta = 1.5;
figure(2) 
plot(time,coh_cos_the(:,2),'k');
hold on
plot(time,work_cos_the(:,2),'--r');
hold on
plot(time,heat_cos_the(:,2),':b');
hold off
xlabel('$\omega t$', 'FontSize', 50, 'Interpreter','latex');
set(gca,'FontSize', 38, 'TickLabelInterpreter','latex')
set(gca,'TickLabelInterpreter','latex')

% Thermal state - beta = 1.5;
figure(3) 
plot(time,entroS_cos_the(:,2),'k');
hold on
plot(time,entroF_cos_the(:,2),'--r');
hold on
plot(time,entroC_cos_the(:,2),':b');
hold off
xlabel('$\omega t$', 'FontSize', 50, 'Interpreter','latex');
set(gca,'FontSize', 38, 'TickLabelInterpreter','latex')
set(gca,'TickLabelInterpreter','latex')

