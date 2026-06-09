% Stephan-Boltzmann law

clear all
close all

sigma = 5.67e-8; % [Jm^-2 s^-1 K^4 ]
T = linspace(0,1000);
S = sigma*T.^4;

rng = find(T>=300 & T<=350);

figure; 
hold on
h(1) = plot(T,S,'r','LineWidth',2);
h(2) = plot(T(rng),S(rng),'k','LineWidth',3);
xlabel('T [K]','Interpreter','latex')
ylabel('S $[Wm^{-2}]$','Interpreter','latex')
grid on
box on
set(gcf,'color','w')
set(gca,'Yscale','lin')
set(gca,'Xscale','lin')
set(gcf,'Position',[360  340 350 300])

figure; 
hold on
h(1) = plot(T,S,'r','LineWidth',2);
h(2) = plot(T(rng),S(rng),'k','LineWidth',2);
xlabel('T [K]','Interpreter','latex')
ylabel('S $[Wm^{-2}]$','Interpreter','latex')
grid on
box on
set(gcf,'color','w')
set(gca,'Yscale','lin')
set(gca,'Xscale','lin')
set(gcf,'Position',[360  340 350 300])
xlim([300,350])


figure; 
hold on
h(1) = plot(S,T,'r','LineWidth',2);
h(2) = plot(S(rng),T(rng),'k','LineWidth',2);
ylabel('T [K]','Interpreter','latex')
xlabel('S $[Wm^{-2}]$','Interpreter','latex')
grid on
box on
set(gcf,'color','w')
set(gca,'Yscale','lin')
set(gca,'Xscale','lin')
set(gcf,'Position',[360  340 350 300])
% ylim([300,350])



